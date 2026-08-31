"""Europe PMC + PubMed E-utilities combined literature collector.

Europe PMC (EMBL-EBI, free API):
  - PubMed / PMC / 学術リポジトリをまとめて全文検索
  - abstractText を直接返すため efetch 不要
  - フィールド検索: GENE_PROTEIN_NAME / TITLE / ABSTRACT / AUTH_ORCID など

PubMed E-utilities (NCBI, public domain):
  - MeSH ターム × 遺伝子シノニムの多段 tier クエリ
  - Europe PMC で取れなかった論文を補完

スコアリング（ Europe PMC + PubMed 共通）:
  4: 公式シンボル × 疾患名 → タイトル両方
  3: 公式シンボル × 疾患名 → タイトル+アブスト
  2: シノニム × 疾患名 → タイトル
  1: シノニム × 疾患名 → アブスト
  0: 片方のみ（除外）

臨床研究を優先し（is_clinical 降順）、その中で関連度 → 年 の順に並べる。
"""
import time
import re
import requests
from typing import Optional

EPMC_API   = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CLINICAL_PUBTYPES = {
    "clinical trial", "clinical trial, phase i", "clinical trial, phase ii",
    "clinical trial, phase iii", "clinical trial, phase iv",
    "controlled clinical trial", "randomized controlled trial",
    "pragmatic clinical trial", "adaptive clinical trial",
    "observational study", "multicenter study", "case reports",
    "journal article",   # EPMC は "Journal Article" が基本 pubtype
}

_CLINICAL_KEYWORDS = {
    "clinical trial", "randomized", "randomised", "placebo", "patients",
    "cohort", "case-control", "observational", "phase i", "phase ii",
    "phase iii", "phase iv",
}


def _is_clinical(pub_types: list[str]) -> bool:
    lower = {(pt or "").lower() for pt in (pub_types or [])}
    return bool(lower & {
        "clinical trial", "clinical trial, phase i", "clinical trial, phase ii",
        "clinical trial, phase iii", "clinical trial, phase iv",
        "controlled clinical trial", "randomized controlled trial",
        "pragmatic clinical trial", "adaptive clinical trial",
        "observational study", "multicenter study", "case reports",
    })


def _is_clinical_by_abstract(abstract: str) -> bool:
    """アブストラクトのキーワードから臨床研究かを補完判定。"""
    al = (abstract or "").lower()
    return any(kw in al for kw in _CLINICAL_KEYWORDS)


# ─────────────────────────────────────────────────────────────────────────────
# シノニム取得（pubmed.py から流用）
# ─────────────────────────────────────────────────────────────────────────────

def _get_gene_synonyms(gene_symbol: str) -> list[str]:
    """NCBI Gene から遺伝子シノニムを取得。公式シンボルを先頭に返す。"""
    synonyms = [gene_symbol]
    seen = {gene_symbol.upper()}

    def _add(name: str):
        name = name.strip()
        if name and name.upper() not in seen and len(name) >= 2:
            seen.add(name.upper())
            synonyms.append(name)

    try:
        r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params={
            "db": "gene",
            "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
            "retmax": 1, "retmode": "json",
        }, timeout=10)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if ids:
            time.sleep(0.2)
            r2 = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params={
                "db": "gene", "id": ids[0],
                "rettype": "gene_table", "retmode": "text",
            }, timeout=10)
            r2.raise_for_status()
            for line in r2.text.splitlines():
                if line.startswith("Also known as"):
                    for s in re.split(r"[;,]", line.replace("Also known as", "").strip()):
                        _add(s)
                    break
    except Exception:
        pass
    return synonyms


def _get_mesh_heading(disease_name: str) -> str:
    """疾患名 → MeSH 見出し語（PubMed クエリ用）。失敗時は disease_name を返す。"""
    try:
        r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params={
            "db": "mesh", "term": disease_name, "retmax": 1, "retmode": "json",
        }, timeout=10)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return disease_name
        time.sleep(0.2)
        r2 = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params={
            "db": "mesh", "id": ids[0],
        }, timeout=10)
        r2.raise_for_status()
        for line in r2.text.splitlines():
            m = re.match(r"^\d+:\s+(.+)$", line.strip())
            if m:
                return m.group(1).strip()
    except Exception:
        pass
    return disease_name


# ─────────────────────────────────────────────────────────────────────────────
# Europe PMC 検索
# ─────────────────────────────────────────────────────────────────────────────

def _epmc_search(query: str, page_size: int = 100, max_pages: int = 3) -> list[dict]:
    """Europe PMC REST API で検索し、論文メタデータのリストを返す。"""
    articles: list[dict] = []
    cursor_mark = "*"
    for _ in range(max_pages):
        try:
            params = {
                "query":       query,
                "resultType":  "core",
                "format":      "json",
                "pageSize":    page_size,
                "cursorMark":  cursor_mark,
                "sort":        "RELEVANCE",
            }
            r = requests.get(EPMC_API, params=params, timeout=25)
            r.raise_for_status()
            data = r.json()
            results = (data.get("resultList") or {}).get("result") or []
            if not results:
                break
            articles.extend(results)
            # 次ページ
            next_cursor = data.get("nextCursorMark")
            if not next_cursor or next_cursor == cursor_mark:
                break
            cursor_mark = next_cursor
        except Exception as e:
            print(f"    [EuropePMC] search error: {e}")
            break
        time.sleep(0.2)
    return articles


def _epmc_to_paper(item: dict) -> dict:
    """Europe PMC の結果アイテム → 統一フォーマットの paper dict。"""
    pmid = str(item.get("pmid") or item.get("id") or "").strip()
    pub_types_raw = (item.get("pubTypeList") or {}).get("pubType") or []
    if isinstance(pub_types_raw, str):
        pub_types_raw = [pub_types_raw]
    authors_raw = (item.get("authorList") or {}).get("author") or []
    if isinstance(authors_raw, list):
        author_names = [
            a.get("fullName") or f"{a.get('lastName', '')} {a.get('initials', '')}".strip()
            for a in authors_raw[:3]
        ]
    else:
        author_names = []

    abstract = (item.get("abstractText") or "").strip()
    abstract = re.sub(r"<[^>]+>", " ", abstract).strip()  # HTML タグ除去

    is_clin = _is_clinical(pub_types_raw) or _is_clinical_by_abstract(abstract)

    return {
        "pmid":            pmid,
        "title":           item.get("title", "").rstrip(".").strip(),
        "journal":         item.get("journalTitle", "") or item.get("journal", {}).get("title", ""),
        "year":            str(item.get("pubYear") or item.get("firstPublicationDate", ""))[:4],
        "authors":         author_names,
        "abstract":        abstract,
        "relevance_score": 0,
        "match_type":      "",
        "pub_types":       pub_types_raw,
        "is_clinical":     is_clin,
        "_source":         "epmc",
    }


# ─────────────────────────────────────────────────────────────────────────────
# PubMed E-utilities 補完検索
# ─────────────────────────────────────────────────────────────────────────────

def _pubmed_supplement(
    gene: str,
    gene_syns: list[str],
    disease: str,
    mesh_heading: str,
    exclude_pmids: set[str],
    max_ids: int = 300,
) -> list[str]:
    """Europe PMC で取れなかった PMID を PubMed から補完して返す。"""

    def _q_term(term: str, field: str = "Title/Abstract") -> str:
        words = term.split()
        if len(words) <= 1:
            return f'"{term}"[{field}]'
        return "(" + " AND ".join(f'"{w}"[{field}]' for w in words) + ")"

    def _q_gene_official():
        return _q_term(gene)

    def _q_gene_syns(max_syn: int = 8):
        terms = [_q_term(s) for s in gene_syns[:max_syn + 1]]
        return "(" + " OR ".join(terms) + ")"

    def _q_disease_mesh():
        heading = mesh_heading or disease
        return f'"{heading}"[MeSH Terms]'

    def _q_disease_text(max_syn: int = 3):
        return _q_term(disease)

    queries = [
        f'"{gene}"[Gene/Protein Name] AND {_q_disease_mesh()}',
        f'{_q_gene_official()} AND {_q_disease_mesh()}',
        f'{_q_gene_official()} AND {_q_disease_text()}',
        f'{_q_gene_syns()} AND {_q_disease_mesh()}',
    ]

    seen: set[str] = set(exclude_pmids)
    new_pmids: list[str] = []
    for q in queries:
        try:
            r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params={
                "db": "pubmed", "term": q,
                "retmax": max_ids, "retmode": "json", "sort": "pub+date",
            }, timeout=15)
            r.raise_for_status()
            for pid in r.json().get("esearchresult", {}).get("idlist", []):
                if pid not in seen:
                    seen.add(pid)
                    new_pmids.append(pid)
        except Exception:
            pass
        time.sleep(0.2)
        if len(new_pmids) >= max_ids:
            break
    return new_pmids


def _pubmed_fetch_papers(pmids: list[str]) -> list[dict]:
    """PMID リストから PubMed esummary + efetch で論文データを返す。"""
    if not pmids:
        return []

    # esummary で書誌情報
    try:
        r = requests.post(f"{EUTILS_BASE}/esummary.fcgi", data={
            "db": "pubmed", "id": ",".join(pmids), "retmode": "json",
        }, timeout=20)
        r.raise_for_status()
        summaries = r.json().get("result", {})
    except Exception:
        summaries = {}

    # efetch でアブストラクト
    time.sleep(0.3)
    abstract_map: dict[str, str] = {}
    try:
        r = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params={
            "db": "pubmed", "id": ",".join(pmids),
            "rettype": "abstract", "retmode": "xml",
        }, timeout=30)
        r.raise_for_status()
        abstract_map = {
            pmid: re.sub(r"<[^>]+>", " ", text).strip()
            for pmid, text in re.findall(
                r"<MedlineCitation[^>]*>.*?<PMID[^>]*>(\d+)</PMID>"
                r".*?<AbstractText[^>]*>(.*?)</AbstractText>"
                r".*?</MedlineCitation>",
                r.text, re.DOTALL,
            )
        }
    except Exception:
        pass

    papers = []
    for pmid in pmids:
        s = summaries.get(pmid) or {}
        pub_types = list(s.get("pubtype") or [])
        abstract  = abstract_map.get(pmid, "")
        is_clin   = _is_clinical(pub_types) or _is_clinical_by_abstract(abstract)
        papers.append({
            "pmid":            pmid,
            "title":           s.get("title", ""),
            "journal":         s.get("fulljournalname", ""),
            "year":            (s.get("pubdate", "") or "")[:4],
            "authors":         [a.get("name", "") for a in (s.get("authors") or [])[:3]],
            "abstract":        abstract,
            "relevance_score": 0,
            "match_type":      "",
            "pub_types":       pub_types,
            "is_clinical":     is_clin,
            "_source":         "pubmed",
        })
    return papers


# ─────────────────────────────────────────────────────────────────────────────
# スコアリング
# ─────────────────────────────────────────────────────────────────────────────

def _score_papers(
    papers: list[dict],
    gene: str,
    gene_syns: list[str],
    disease: str,
):
    """タイトル + アブストラクトのテキストマッチでスコアを付与（in-place）。"""
    official_l   = gene.lower()
    gene_syns_l  = [s.lower() for s in gene_syns[1:]]
    disease_l    = disease.lower()
    disease_words = [w for w in disease_l.split() if len(w) > 3]

    def _d_match(text: str) -> bool:
        t = text.lower()
        return disease_l in t or all(w in t for w in disease_words)

    for p in papers:
        tl = p["title"].lower()
        al = (p["abstract"] or "").lower()

        off_t  = official_l in tl
        off_a  = official_l in al
        dis_t  = _d_match(p["title"])
        dis_a  = _d_match(p["abstract"] or "")
        syn_t  = any(s in tl for s in gene_syns_l)
        syn_a  = any(s in al for s in gene_syns_l)

        if off_t and dis_t:
            score, match = 4, "公式シンボル×タイトル"
        elif (off_t or off_a) and (dis_t or dis_a):
            score, match = 3, "公式シンボル×アブスト"
        elif syn_t and dis_t:
            score, match = 2, "シノニム×タイトル"
        elif (syn_t or syn_a) and (dis_t or dis_a):
            score, match = 1, "シノニム×アブスト"
        else:
            score, match = 0, "内容参照"

        p["relevance_score"] = score
        p["match_type"]      = match


# ─────────────────────────────────────────────────────────────────────────────
# メイン関数（search_pubmed と同じシグネチャ・出力フォーマット）
# ─────────────────────────────────────────────────────────────────────────────

def search_literature(
    gene: str,
    disease: str,
    max_results: int = 100,
    disease_efo_id: str = None,
) -> list[dict]:
    """Europe PMC + PubMed E-utilities から論文を収集してスコア順に返す。

    出力フォーマット（search_pubmed と同一）:
      pmid, title, journal, year, authors, abstract,
      relevance_score, match_type, pub_types, is_clinical
    """
    # ── シノニム + MeSH 見出し語 ──────────────────────────────────────────────
    gene_syns   = _get_gene_synonyms(gene)
    mesh_heading = _get_mesh_heading(disease)

    n_gsyn = len(gene_syns) - 1
    print(f"    [文献] 遺伝子シノニム ({len(gene_syns)}件): "
          f"{', '.join(gene_syns[:6])}{'...' if n_gsyn > 5 else ''}")
    print(f"    [文献] MeSH: 「{mesh_heading}」")

    # ── Europe PMC 検索 ────────────────────────────────────────────────────
    # フィールド指定クエリ: 遺伝子シンボル × 疾患名でタイトル+アブスト絞り込み
    gene_q    = " OR ".join(f'"{s}"' for s in gene_syns[:5])
    disease_q = f'"{disease}"'
    if mesh_heading and mesh_heading.lower() != disease.lower():
        disease_q = f'("{disease}" OR "{mesh_heading}")'

    epmc_query = (
        f'(TITLE:({gene_q}) OR ABSTRACT:({gene_q})) AND '
        f'(TITLE:{disease_q} OR ABSTRACT:{disease_q})'
    )
    raw_epmc = _epmc_search(epmc_query, page_size=100, max_pages=3)
    print(f"    [文献] Europe PMC ヒット: {len(raw_epmc)} 件")

    # Europe PMC → 統一フォーマット変換（PMID ありのみ）
    epmc_papers: list[dict] = []
    seen_pmids:  set[str]   = set()
    for item in raw_epmc:
        p = _epmc_to_paper(item)
        if p["pmid"] and p["pmid"] not in seen_pmids:
            seen_pmids.add(p["pmid"])
            epmc_papers.append(p)

    # ── PubMed 補完 ────────────────────────────────────────────────────────
    supplement_pmids = _pubmed_supplement(
        gene, gene_syns, disease, mesh_heading,
        exclude_pmids=seen_pmids,
        max_ids=max_results * 3,
    )
    time.sleep(0.3)
    pubmed_papers = _pubmed_fetch_papers(supplement_pmids[:300])
    print(f"    [文献] PubMed 補完: {len(pubmed_papers)} 件")

    # ── マージ ─────────────────────────────────────────────────────────────
    all_papers = epmc_papers + pubmed_papers

    # ── スコアリング ───────────────────────────────────────────────────────
    _score_papers(all_papers, gene, gene_syns, disease)

    # score == 0 は除外
    all_papers = [p for p in all_papers if p["relevance_score"] > 0]

    # Europe PMC は relevance ソート済みなので同スコア内で EPMC を優先
    def _sort_key(p):
        src_priority = 0 if p.get("_source") == "epmc" else 1
        return (p["is_clinical"], p["relevance_score"], p["year"], -src_priority)

    all_papers.sort(key=_sort_key, reverse=True)

    # _source フィールドは内部用なので除去
    for p in all_papers:
        p.pop("_source", None)

    return all_papers[:max_results]
