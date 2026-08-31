"""PubTator3 biomedical entity-annotated literature collector (public domain).

PubTator3 は PubMed 論文に遺伝子・疾患・化合物・変異などのエンティティを
自動アノテーションしたデータベース。

検索戦略（優先順）:
  Tier 1: NCBI Gene ID × MeSH Disease ID でエンティティ共起論文を取得
  Tier 2: 遺伝子シンボル × 疾患名のテキスト検索（ID解決失敗・補完用）
  Tier 3: PubMed E-utilities フォールバック（PubTator3 が 0 件の場合）

スコアリング（BioC JSON アノテーションベース）:
  4: 遺伝子エンティティ × 疾患エンティティ → タイトル両方に存在
  3: 遺伝子エンティティ × 疾患エンティティ → タイトル+アブスト
  2: 遺伝子テキスト × 疾患テキスト → タイトル
  1: 遺伝子テキスト × 疾患テキスト → アブスト
  0: 片方のみ（除外）
"""
import time
import re
import requests
from typing import Optional

PUBTATOR_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
EUTILS_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CLINICAL_PUBTYPES = {
    "clinical trial", "clinical trial, phase i", "clinical trial, phase ii",
    "clinical trial, phase iii", "clinical trial, phase iv",
    "controlled clinical trial", "randomized controlled trial",
    "pragmatic clinical trial", "adaptive clinical trial",
    "observational study", "multicenter study", "case reports",
}


def _is_clinical(pub_types: list[str]) -> bool:
    return any((pt or "").lower() in CLINICAL_PUBTYPES for pt in (pub_types or []))


# ─────────────────────────────────────────────────────────────────────────────
# ID 解決
# ─────────────────────────────────────────────────────────────────────────────

def _resolve_ncbi_gene_id(gene_symbol: str) -> Optional[str]:
    """遺伝子シンボル → NCBI Gene ID (ヒト)。失敗時は None。"""
    try:
        r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params={
            "db": "gene",
            "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
            "retmax": 1, "retmode": "json",
        }, timeout=10)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        return ids[0] if ids else None
    except Exception:
        return None


def _resolve_mesh_id(disease_name: str) -> Optional[str]:
    """疾患名 → MeSH UI (例: D000544)。失敗時は None。"""
    try:
        r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params={
            "db": "mesh", "term": disease_name,
            "retmax": 1, "retmode": "json",
        }, timeout=10)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None
        time.sleep(0.2)
        r2 = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params={
            "db": "mesh", "id": ids[0],
        }, timeout=10)
        r2.raise_for_status()
        # "MeSH Unique ID:" 行から D-number を抽出
        m = re.search(r"MeSH Unique ID:\s*(D\d+)", r2.text)
        return m.group(1) if m else None
    except Exception:
        return None


# ─────────────────────────────────────────────────────────────────────────────
# PubTator3 検索
# ─────────────────────────────────────────────────────────────────────────────

def _pubtator_search(query: str, max_pages: int = 3) -> list[str]:
    """PubTator3 search API で PMID リストを返す（重複除去済み）。

    PubTator3 API レスポンス形式:
      {"count": N, "total": M, "results": [{pmid: int, ...}, ...]}
    """
    pmids: list[str] = []
    seen: set[str] = set()
    for page in range(1, max_pages + 1):
        try:
            r = requests.get(f"{PUBTATOR_BASE}/search/", params={
                "text": query, "page": page,
            }, timeout=20)
            r.raise_for_status()
            data = r.json()
            # API は "results" キーを使用（旧APIの "hits" も念のため参照）
            hits = data.get("results") or data.get("hits") or []
            if not hits:
                break
            for h in hits:
                # pmid は integer または string で返る
                pmid = str(h.get("pmid") or h.get("_id") or "").strip()
                if pmid and pmid != "0" and pmid not in seen:
                    seen.add(pmid)
                    pmids.append(pmid)
            # 次ページがなければ終了
            total = data.get("total", 0)
            if len(pmids) >= total:
                break
        except Exception as e:
            print(f"    [PubTator] search error (query={query!r}): {e}")
            break
        time.sleep(0.3)
    return pmids


# ─────────────────────────────────────────────────────────────────────────────
# BioC JSON エクスポート（アノテーション付きアブストラクト）
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_biocjson(pmids: list[str]) -> dict[str, dict]:
    """PMID リストの BioC JSON を取得し {pmid: doc} を返す。"""
    if not pmids:
        return {}
    result: dict[str, dict] = {}
    batch_size = 100
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        try:
            r = requests.get(f"{PUBTATOR_BASE}/publications/export/biocjson", params={
                "pmids": ",".join(batch),
            }, timeout=30)
            r.raise_for_status()
            raw = r.json()
            # レスポンスは配列 or {"PubTator3": [...]} のどちらかで返ることがある
            if isinstance(raw, list):
                docs = raw
            elif isinstance(raw, dict):
                docs = raw.get("PubTator3") or raw.get("results") or [raw]
            else:
                docs = []
            for doc in docs:
                # id フィールドは "pmid" または "_id" または "id" で来ることがある
                pmid = str(
                    doc.get("pmid") or doc.get("_id") or doc.get("id") or ""
                ).strip()
                if pmid:
                    result[pmid] = doc
        except Exception as e:
            print(f"    [PubTator] biocjson error: {e}")
        time.sleep(0.3)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# NCBI esummary（pub_types / 書誌情報）
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_esummary(pmids: list[str]) -> dict[str, dict]:
    """PMID → esummary (pub_types, title, journal, year, authors)。"""
    if not pmids:
        return {}
    try:
        r = requests.post(f"{EUTILS_BASE}/esummary.fcgi", data={
            "db": "pubmed", "id": ",".join(pmids), "retmode": "json",
        }, timeout=20)
        r.raise_for_status()
        return r.json().get("result", {})
    except Exception:
        return {}


# ─────────────────────────────────────────────────────────────────────────────
# PubMed E-utilities フォールバック検索（PubTator3 が 0 件の場合に使用）
# ─────────────────────────────────────────────────────────────────────────────

def _pubmed_fallback_search(gene: str, disease: str, max_results: int) -> list[str]:
    """PubTator3 が 0 件の場合に PubMed esearch を使って PMID を返す。"""
    queries = [
        f'"{gene}"[Title/Abstract] AND "{disease}"[Title/Abstract]',
        f'"{gene}"[Gene/Protein Name] AND "{disease}"[MeSH Terms]',
        f'"{gene}"[Title/Abstract] AND "{disease}"[MeSH Terms]',
    ]
    seen: set[str] = set()
    pmids: list[str] = []
    for q in queries:
        try:
            r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params={
                "db": "pubmed", "term": q,
                "retmax": max_results * 2,
                "retmode": "json", "sort": "pub+date",
            }, timeout=15)
            r.raise_for_status()
            for pid in r.json().get("esearchresult", {}).get("idlist", []):
                if pid not in seen:
                    seen.add(pid)
                    pmids.append(pid)
            if len(pmids) >= max_results * 3:
                break
        except Exception:
            pass
        time.sleep(0.2)
    return pmids


# ─────────────────────────────────────────────────────────────────────────────
# アブストラクト取得（esummary で取れない場合の補完）
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_abstracts_efetch(pmids: list[str]) -> dict[str, str]:
    """efetch XML から {pmid: abstract} を返す。"""
    if not pmids:
        return {}
    try:
        r = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params={
            "db": "pubmed", "id": ",".join(pmids),
            "rettype": "abstract", "retmode": "xml",
        }, timeout=30)
        r.raise_for_status()
        return {
            pmid: re.sub(r"<[^>]+>", " ", text).strip()
            for pmid, text in re.findall(
                r"<MedlineCitation[^>]*>.*?<PMID[^>]*>(\d+)</PMID>"
                r".*?<AbstractText[^>]*>(.*?)</AbstractText>"
                r".*?</MedlineCitation>",
                r.text, re.DOTALL,
            )
        }
    except Exception:
        return {}


# ─────────────────────────────────────────────────────────────────────────────
# スコアリング（アノテーションベース + テキストフォールバック）
# ─────────────────────────────────────────────────────────────────────────────

def _score_doc(
    doc: dict,
    gene_symbol: str,
    ncbi_gene_id: Optional[str],
    disease_name: str,
    mesh_id: Optional[str],
    abstract_text: str = "",
    title_text: str = "",
) -> tuple[int, str]:
    """BioC JSON doc から関連度スコアと match_type を返す。
    doc が空の場合は title_text / abstract_text のテキストマッチを使う。
    """
    gene_upper    = gene_symbol.upper()
    disease_lower = disease_name.lower()
    disease_words = [w for w in disease_lower.split() if len(w) > 3]

    # BioC JSON が取れている場合はアノテーションベース
    if doc:
        title_annots:    list[dict] = []
        abstract_annots: list[dict] = []
        doc_title_text = ""
        doc_abstract_text = ""

        for passage in (doc.get("passages") or []):
            section = (passage.get("infons", {}).get("section_type") or
                       passage.get("infons", {}).get("type") or "").upper()
            annots = passage.get("annotations") or []
            text   = passage.get("text", "")
            if section in ("TITLE", "FRONT"):
                title_annots.extend(annots)
                doc_title_text += " " + text
            else:
                abstract_annots.extend(annots)
                doc_abstract_text += " " + text

        # アノテーション内のターゲット遺伝子・疾患を確認
        def _has_gene_annot(annots: list[dict]) -> bool:
            for a in annots:
                t = (a.get("infons", {}).get("type") or "").lower()
                if t != "gene":
                    continue
                name = (a.get("text") or "").upper()
                ids_raw = str(a.get("infons", {}).get("identifier") or "")
                ids = [x.strip() for x in ids_raw.split(",") if x.strip()]
                if name == gene_upper:
                    return True
                if ncbi_gene_id and ncbi_gene_id in ids:
                    return True
            return False

        def _has_disease_annot(annots: list[dict]) -> bool:
            for a in annots:
                t = (a.get("infons", {}).get("type") or "").lower()
                if t != "disease":
                    continue
                name = (a.get("text") or "").lower()
                ids_raw = str(a.get("infons", {}).get("identifier") or "")
                ids = [x.strip() for x in ids_raw.split(",") if x.strip()]
                if mesh_id and (f"MESH:{mesh_id}" in ids or mesh_id in ids):
                    return True
                if any(w in name for w in disease_words):
                    return True
            return False

        g_t = _has_gene_annot(title_annots)
        g_a = _has_gene_annot(abstract_annots)
        d_t = _has_disease_annot(title_annots)
        d_a = _has_disease_annot(abstract_annots)

        if g_t and d_t:
            return 4, "エンティティ×タイトル"
        if (g_t or g_a) and (d_t or d_a):
            return 3, "エンティティ×アブスト"

        # テキストフォールバック（doc 内テキスト使用）
        title_text    = title_text    or doc_title_text
        abstract_text = abstract_text or doc_abstract_text

    # テキストベーススコアリング
    tl = title_text.upper()
    al = abstract_text.upper()
    dl = disease_lower

    def _d_in(text: str) -> bool:
        t = text.lower()
        return dl in t or all(w in t for w in disease_words)

    g_in_title = gene_upper in tl
    g_in_abs   = gene_upper in al
    d_in_title = _d_in(title_text)
    d_in_abs   = _d_in(abstract_text)

    if g_in_title and d_in_title:
        return 2, "テキスト×タイトル"
    if (g_in_title or g_in_abs) and (d_in_title or d_in_abs):
        return 1, "テキスト×アブスト"
    return 0, "内容参照"


# ─────────────────────────────────────────────────────────────────────────────
# メイン関数（pubmed.search_pubmed と同じシグネチャ・出力フォーマット）
# ─────────────────────────────────────────────────────────────────────────────

def search_pubtator(
    gene: str,
    disease: str,
    max_results: int = 100,
    disease_efo_id: str = None,
) -> list[dict]:
    """Return top PubTator3-annotated abstracts scored by entity co-occurrence.

    pubmed.search_pubmed と同じ出力フォーマット:
      pmid, title, journal, year, authors, abstract,
      relevance_score, match_type, pub_types, is_clinical
    """
    # ── ID 解決 ──────────────────────────────────────────────────────────────
    ncbi_gene_id = _resolve_ncbi_gene_id(gene)
    mesh_id      = _resolve_mesh_id(disease)
    print(f"    PubTator: NCBI Gene ID={ncbi_gene_id}, MeSH ID={mesh_id}")

    # ── PubTator3 検索 ───────────────────────────────────────────────────────
    all_pmids: list[str] = []
    seen: set[str] = set()

    def _add(pmids):
        for p in pmids:
            if p not in seen:
                seen.add(p)
                all_pmids.append(p)

    # Tier 1: NCBI Gene ID × MeSH Disease ID（最高精度）
    # PubTator3 概念検索: "@{ncbi_gene_id}" と "@MESH:{mesh_id}" の形式
    if ncbi_gene_id and mesh_id:
        _add(_pubtator_search(f"@{ncbi_gene_id} @MESH:{mesh_id}", max_pages=5))

    # Tier 2: Gene ID のみ（mesh_id が取れない場合）
    if ncbi_gene_id and len(all_pmids) < max_results:
        _add(_pubtator_search(f"@{ncbi_gene_id} {disease}", max_pages=3))

    # Tier 3: テキスト検索（補完 or ID解決失敗時）
    if len(all_pmids) < max_results:
        _add(_pubtator_search(f"{gene} {disease}", max_pages=3))

    # Tier 4: PubMed E-utilities フォールバック（PubTator3 が 0 件の場合）
    if not all_pmids:
        print(f"    PubTator: 0 hits, falling back to PubMed E-utilities...")
        _add(_pubmed_fallback_search(gene, disease, max_results))

    if not all_pmids:
        return []

    print(f"    PubTator: {len(all_pmids)} unique PMIDs collected")

    # ── BioC JSON エクスポート（アノテーション付き） ───────────────────────
    fetch_pmids = all_pmids[:max_results * 3]
    bioc_docs   = _fetch_biocjson(fetch_pmids)

    # ── esummary（書誌情報 + pub_types） ─────────────────────────────────────
    time.sleep(0.3)
    summary = _fetch_esummary(fetch_pmids)

    # BioC JSON でアブストラクトが取れなかった PMID を efetch で補完
    missing_abs = [p for p in fetch_pmids if p not in bioc_docs or not any(
        (pa.get("infons", {}).get("section_type") or "").upper() not in ("TITLE", "FRONT")
        for pa in (bioc_docs.get(p) or {}).get("passages", [])
    )]
    efetch_abstracts: dict[str, str] = {}
    if missing_abs:
        time.sleep(0.3)
        efetch_abstracts = _fetch_abstracts_efetch(missing_abs[:200])

    # ── 論文データ構築 + スコアリング ─────────────────────────────────────────
    papers: list[dict] = []
    for pmid in fetch_pmids:
        s   = summary.get(pmid) or {}
        doc = bioc_docs.get(pmid) or {}

        # タイトル（esummary 優先、なければ BioC JSON）
        title = s.get("title", "")
        if not title and doc.get("passages"):
            for p in doc["passages"]:
                if (p.get("infons", {}).get("section_type") or "").upper() in ("TITLE", "FRONT"):
                    title = p.get("text", "")
                    break

        # アブストラクト（BioC JSON passages から結合）
        abstract_parts = []
        for passage in (doc.get("passages") or []):
            section = (passage.get("infons", {}).get("section_type") or
                       passage.get("infons", {}).get("type") or "").upper()
            if section not in ("TITLE", "FRONT"):
                abstract_parts.append(passage.get("text", ""))
        abstract = " ".join(abstract_parts).strip()
        if not abstract:
            abstract = efetch_abstracts.get(pmid, "")

        pub_types = list(s.get("pubtype") or [])
        score, match = _score_doc(
            doc, gene, ncbi_gene_id, disease, mesh_id,
            abstract_text=abstract, title_text=title,
        )

        papers.append({
            "pmid":            pmid,
            "title":           title,
            "journal":         s.get("fulljournalname", ""),
            "year":            (s.get("pubdate", "") or "")[:4],
            "authors":         [a.get("name", "") for a in (s.get("authors") or [])[:3]],
            "abstract":        abstract,
            "relevance_score": score,
            "match_type":      match,
            "pub_types":       pub_types,
            "is_clinical":     _is_clinical(pub_types),
        })

    # score == 0 は除外
    papers = [p for p in papers if p["relevance_score"] > 0]

    # 臨床研究優先 → 関連度 → 年の降順
    papers.sort(key=lambda p: (p["is_clinical"], p["relevance_score"], p["year"]), reverse=True)

    return papers[:max_results]


def fetch_abstract(pmid: str) -> Optional[str]:
    """単一 PMID のアブストラクトを返す（後方互換）。"""
    docs = _fetch_biocjson([pmid])
    doc  = docs.get(pmid)
    if doc:
        parts = []
        for passage in (doc.get("passages") or []):
            section = (passage.get("infons", {}).get("section_type") or "").upper()
            if section not in ("TITLE", "FRONT"):
                parts.append(passage.get("text", ""))
        text = " ".join(parts).strip()
        if text:
            return text
    # フォールバック: efetch
    result = _fetch_abstracts_efetch([pmid])
    return result.get(pmid)
