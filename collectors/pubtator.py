"""PubTator3 biomedical entity-annotated literature collector (public domain).

PubTator3 は PubMed 論文に遺伝子・疾患・化合物・変異などのエンティティを
自動アノテーションしたデータベース。

検索戦略:
  Tier 1: NCBI Gene ID × MeSH Disease ID でエンティティ共起論文を取得（最高精度）
  Tier 2: 遺伝子シンボル × 疾患名テキストで追加取得（ID解決失敗時のフォールバック）

スコアリング（BioC JSON アノテーションベース）:
  4: 遺伝子エンティティ × 疾患エンティティ → タイトル両方
  3: 遺伝子エンティティ × 疾患エンティティ → タイトル+アブストラクト
  2: 遺伝子テキスト × 疾患テキスト → タイトル
  1: 遺伝子テキスト × 疾患テキスト → アブストラクト
  0: どちらか一方のみ（除外対象）
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
# ID解決
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
    """PubTator3 search API で PMID リストを返す（重複除去済み）。"""
    pmids: list[str] = []
    seen: set[str] = set()
    for page in range(1, max_pages + 1):
        try:
            r = requests.get(f"{PUBTATOR_BASE}/search/", params={
                "text": query, "page": page,
            }, timeout=20)
            r.raise_for_status()
            data = r.json()
            hits = data.get("hits") or []
            if not hits:
                break
            for h in hits:
                pmid = str(h.get("pmid", ""))
                if pmid and pmid not in seen:
                    seen.add(pmid)
                    pmids.append(pmid)
        except Exception:
            break
        time.sleep(0.2)
    return pmids


# ─────────────────────────────────────────────────────────────────────────────
# BioC JSON エクスポート（アノテーション付きアブストラクト）
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_biocjson(pmids: list[str]) -> dict[str, dict]:
    """PMID リストの BioC JSON を取得し {pmid: doc} を返す。"""
    if not pmids:
        return {}
    result: dict[str, dict] = {}
    # 一度に 100 件まで
    batch_size = 100
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        try:
            r = requests.get(f"{PUBTATOR_BASE}/publications/export/biocjson", params={
                "pmids": ",".join(batch),
            }, timeout=30)
            r.raise_for_status()
            data = r.json()
            docs = data if isinstance(data, list) else data.get("PubTator3", [data])
            for doc in docs:
                pmid = str(doc.get("pmid", "") or doc.get("id", ""))
                if pmid:
                    result[pmid] = doc
        except Exception:
            pass
        time.sleep(0.3)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# NCBI esummary（pub_types / 書誌情報 取得）
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_esummary(pmids: list[str]) -> dict[str, dict]:
    """PMID → esummary (pub_types, title, journal, year, authors) を返す。"""
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
# スコアリング（アノテーションベース）
# ─────────────────────────────────────────────────────────────────────────────

def _score_doc(
    doc: dict,
    gene_symbol: str,
    ncbi_gene_id: Optional[str],
    disease_name: str,
    mesh_id: Optional[str],
) -> tuple[int, str]:
    """BioC JSON doc から関連度スコアと match_type を返す。"""
    gene_upper = gene_symbol.upper()
    disease_lower = disease_name.lower()
    disease_words = [w for w in disease_lower.split() if len(w) > 3]

    # パッセージをタイトル / アブストラクトに分類
    title_annots:    list[dict] = []
    abstract_annots: list[dict] = []
    title_text = ""
    abstract_text = ""

    for passage in (doc.get("passages") or []):
        section = (passage.get("infons", {}).get("section_type") or
                   passage.get("infons", {}).get("type") or "").upper()
        annots = passage.get("annotations") or []
        text   = passage.get("text", "")
        if section in ("TITLE", "FRONT"):
            title_annots.extend(annots)
            title_text += " " + text
        else:
            abstract_annots.extend(annots)
            abstract_text += " " + text

    def _has_gene_annot(annots: list[dict]) -> bool:
        for a in annots:
            t = (a.get("infons", {}).get("type") or "").lower()
            if t != "gene":
                continue
            name = (a.get("text") or "").upper()
            ids_raw = a.get("infons", {}).get("identifier") or ""
            ids = [x.strip() for x in str(ids_raw).split(",") if x.strip()]
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
            ids_raw = a.get("infons", {}).get("identifier") or ""
            ids = [x.strip() for x in str(ids_raw).split(",") if x.strip()]
            if mesh_id and f"MESH:{mesh_id}" in ids:
                return True
            if mesh_id and mesh_id in ids:
                return True
            if any(w in name for w in disease_words):
                return True
        return False

    def _text_has_gene(text: str) -> bool:
        return gene_upper in text.upper()

    def _text_has_disease(text: str) -> bool:
        tl = text.lower()
        return disease_lower in tl or all(w in tl for w in disease_words)

    gene_in_title_annot    = _has_gene_annot(title_annots)
    gene_in_abs_annot      = _has_gene_annot(abstract_annots)
    disease_in_title_annot = _has_disease_annot(title_annots)
    disease_in_abs_annot   = _has_disease_annot(abstract_annots)

    gene_in_title_text    = _text_has_gene(title_text)
    gene_in_abs_text      = _text_has_gene(abstract_text)
    disease_in_title_text = _text_has_disease(title_text)
    disease_in_abs_text   = _text_has_disease(abstract_text)

    # アノテーションベース（高精度）
    if gene_in_title_annot and disease_in_title_annot:
        return 4, "エンティティ×タイトル"
    if (gene_in_title_annot or gene_in_abs_annot) and (disease_in_title_annot or disease_in_abs_annot):
        return 3, "エンティティ×アブスト"
    # テキストベース（フォールバック）
    if gene_in_title_text and disease_in_title_text:
        return 2, "テキスト×タイトル"
    if (gene_in_title_text or gene_in_abs_text) and (disease_in_title_text or disease_in_abs_text):
        return 1, "テキスト×アブスト"
    return 0, "内容参照"


# ─────────────────────────────────────────────────────────────────────────────
# メイン関数（pubmed.search_pubmed と同じシグネチャ）
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

    # ── 検索（Tier 1: ID直接指定 → Tier 2: テキスト） ───────────────────────
    all_pmids: list[str] = []
    seen: set[str] = set()

    def _add(pmids):
        for p in pmids:
            if p not in seen:
                seen.add(p)
                all_pmids.append(p)

    if ncbi_gene_id and mesh_id:
        # Tier 1: NCBI Gene ID × MeSH Disease ID（最高精度）
        _add(_pubtator_search(f"@{ncbi_gene_id} @MESH:{mesh_id}", max_pages=5))

    if len(all_pmids) < max_results:
        # Tier 2: シンボル × 疾患名テキスト（補完）
        _add(_pubtator_search(f"{gene} {disease}", max_pages=3))

    if not ncbi_gene_id or not mesh_id:
        # ID解決失敗時の追加フォールバック
        _add(_pubtator_search(f'"{gene}" "{disease}"', max_pages=2))

    if not all_pmids:
        return []

    # ── BioC JSON エクスポート（アノテーション付き） ───────────────────────
    fetch_pmids = all_pmids[:max_results * 3]
    bioc_docs   = _fetch_biocjson(fetch_pmids)

    # ── esummary（書誌情報 + pub_types） ─────────────────────────────────────
    time.sleep(0.3)
    summary = _fetch_esummary(fetch_pmids)

    # ── 論文データ構築 + スコアリング ─────────────────────────────────────────
    papers: list[dict] = []
    for pmid in fetch_pmids:
        s    = summary.get(pmid) or {}
        doc  = bioc_docs.get(pmid) or {}

        # アブストラクトテキスト（BioC JSON passages から結合）
        abstract_parts = []
        for passage in (doc.get("passages") or []):
            section = (passage.get("infons", {}).get("section_type") or
                       passage.get("infons", {}).get("type") or "").upper()
            if section not in ("TITLE", "FRONT"):
                abstract_parts.append(passage.get("text", ""))
        abstract = " ".join(abstract_parts).strip()

        pub_types = list(s.get("pubtype") or [])
        score, match = _score_doc(doc, gene, ncbi_gene_id, disease, mesh_id)

        papers.append({
            "pmid":            pmid,
            "title":           s.get("title", "") or (doc.get("passages") or [{}])[0].get("text", ""),
            "journal":         s.get("fulljournalname", ""),
            "year":            (s.get("pubdate", "") or "")[:4],
            "authors":         [a.get("name", "") for a in (s.get("authors") or [])[:3]],
            "abstract":        abstract,
            "relevance_score": score,
            "match_type":      match,
            "pub_types":       pub_types,
            "is_clinical":     _is_clinical(pub_types),
        })

    # score == 0（間接一致のみ）は除外
    papers = [p for p in papers if p["relevance_score"] > 0]

    # 臨床研究優先、次いで関連度・年の降順
    papers.sort(key=lambda p: (p["is_clinical"], p["relevance_score"], p["year"]), reverse=True)

    return papers[:max_results]


def fetch_abstract(pmid: str) -> Optional[str]:
    """単一 PMID のアブストラクトを返す（後方互換）。"""
    docs = _fetch_biocjson([pmid])
    doc  = docs.get(pmid)
    if not doc:
        return None
    parts = []
    for passage in (doc.get("passages") or []):
        section = (passage.get("infons", {}).get("section_type") or "").upper()
        if section not in ("TITLE", "FRONT"):
            parts.append(passage.get("text", ""))
    return " ".join(parts).strip() or None
