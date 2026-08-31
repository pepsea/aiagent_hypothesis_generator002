"""Europe PMC + PubMed E-utilities combined literature collector.

Europe PMC (EMBL-EBI, free API, Apache 2.0 compatible):
  - PubMed / PMC / 学術リポジトリをまとめて検索
  - abstractText を直接返すため efetch 不要
  - シンプルなフリーテキストクエリで確実にヒット

PubMed E-utilities (NCBI, public domain):
  - MeSH × 遺伝子シノニムの多段 tier クエリで Europe PMC を補完

スコアリング:
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

EPMC_API    = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

_CLINICAL_PUBTYPES = {
    "clinical trial", "clinical trial, phase i", "clinical trial, phase ii",
    "clinical trial, phase iii", "clinical trial, phase iv",
    "controlled clinical trial", "randomized controlled trial",
    "pragmatic clinical trial", "adaptive clinical trial",
    "observational study", "multicenter study", "case reports",
}
_CLINICAL_KEYWORDS = {
    "clinical trial", "randomized", "randomised", "placebo", "patients",
    "cohort", "case-control", "observational", "phase i", "phase ii",
    "phase iii", "phase iv", "double-blind",
}


def _is_clinical(pub_types: list[str], abstract: str = "") -> bool:
    lower = {(pt or "").lower() for pt in (pub_types or [])}
    if lower & _CLINICAL_PUBTYPES:
        return True
    al = (abstract or "").lower()
    return any(kw in al for kw in _CLINICAL_KEYWORDS)


# ─────────────────────────────────────────────────────────────────────────────
# シノニム取得
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
    """疾患名 → MeSH 見出し語。失敗時は disease_name をそのまま返す。"""
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
                "query":      query,
                "resultType": "core",
                "format":     "json",
                "pageSize":   page_size,
                "cursorMark": cursor_mark,
            }
            r = requests.get(EPMC_API, params=params, timeout=25)
            r.raise_for_status()
            data = r.json()
            results = (data.get("resultList") or {}).get("result") or []
            if not results:
                break
            articles.extend(results)
            next_cursor = data.get("nextCursorMark")
            if not next_cursor or next_cursor == cursor_mark:
                break
            cursor_mark = next_cursor
            if len(articles) >= page_size * max_pages:
                break
        except Exception as e:
            print(f"    [EuropePMC] search error (query={repr(query)[:60]}): {e}")
            break
        time.sleep(0.25)
    return articles


def _epmc_to_paper(item: dict) -> dict:
    """Europe PMC アイテム → 統一フォーマット。
    プレプリント (PPR) など PMID がない論文は Europe PMC の id を使う。
    """
    pmid    = str(item.get("pmid")   or "").strip()
    pmcid   = str(item.get("pmcid")  or "").strip()
    src     = str(item.get("source") or "").strip().upper()
    art_id  = str(item.get("id")     or "").strip()

    # 識別子: PMID > PMCID > article_id（PPR1217281 など）
    effective_pmid = pmid or pmcid or art_id

    # リンク URL
    if pmid:
        article_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    elif pmcid:
        article_url = f"https://europepmc.org/article/PMC/{pmcid.replace('PMC','')}"
    elif art_id:
        article_url = f"https://europepmc.org/article/{src}/{art_id}"
    else:
        article_url = ""

    pub_types_raw = (item.get("pubTypeList") or {}).get("pubType") or []
    if isinstance(pub_types_raw, str):
        pub_types_raw = [pub_types_raw]

    authors_raw = (item.get("authorList") or {}).get("author") or []
    if isinstance(authors_raw, list):
        author_names = [
            a.get("fullName") or
            f"{a.get('lastName', '')} {a.get('initials', '')}".strip()
            for a in authors_raw[:3]
        ]
    else:
        author_names = []

    abstract = re.sub(r"<[^>]+>", " ", item.get("abstractText") or "").strip()
    is_clin = _is_clinical(pub_types_raw, abstract)

    return {
        "pmid":            effective_pmid,
        "url":             article_url,
        "title":           (item.get("title") or "").rstrip(".").strip(),
        "journal":         item.get("journalTitle") or "",
        "year":            str(item.get("pubYear") or
                              (item.get("firstPublicationDate") or "")[:4])[:4],
        "authors":         author_names,
        "abstract":        abstract,
        "relevance_score": 0,
        "match_type":      "",
        "pub_types":       pub_types_raw,
        "is_clinical":     is_clin,
        "_source":         "epmc",
    }


# ─────────────────────────────────────────────────────────────────────────────
# PubMed E-utilities 補完
# ─────────────────────────────────────────────────────────────────────────────

def _pubmed_supplement(
    gene: str,
    gene_syns: list[str],
    disease: str,
    mesh_heading: str,
    exclude_pmids: set[str],
    max_ids: int = 300,
) -> list[str]:
    """PubMed esearch で Europe PMC にない PMID を補完して返す。"""

    def _qt(term: str, field: str = "Title/Abstract") -> str:
        words = term.split()
        if len(words) <= 1:
            return f'"{term}"[{field}]'
        return "(" + " AND ".join(f'"{w}"[{field}]' for w in words) + ")"

    mesh_q    = f'"{mesh_heading or disease}"[MeSH Terms]'
    gene_q    = _qt(gene)
    syns_q    = "(" + " OR ".join(_qt(s) for s in gene_syns[:6]) + ")"

    queries = [
        f'"{gene}"[Gene/Protein Name] AND {mesh_q}',
        f'{gene_q} AND {mesh_q}',
        f'{gene_q} AND {_qt(disease)}',
        f'{syns_q} AND {mesh_q}',
        f'{syns_q} AND {_qt(disease)}',
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
    """PMID リストから esummary + efetch で論文データを返す。"""
    if not pmids:
        return []

    try:
        r = requests.post(f"{EUTILS_BASE}/esummary.fcgi", data={
            "db": "pubmed", "id": ",".join(pmids), "retmode": "json",
        }, timeout=20)
        r.raise_for_status()
        summaries = r.json().get("result", {})
    except Exception:
        summaries = {}

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
            "is_clinical":     _is_clinical(pub_types, abstract),
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
    disease_alts: list[str],
):
    """タイトル + アブストラクトのテキストマッチでスコアを付与（in-place）。"""
    official_l    = gene.lower()
    gene_syns_l   = [s.lower() for s in gene_syns[1:]]
    disease_lower_list = [d.lower() for d in disease_alts]
    disease_words = [w for w in disease.lower().split() if len(w) > 3]

    def _d_match(text: str) -> bool:
        t = text.lower()
        return (any(d in t for d in disease_lower_list)
                or all(w in t for w in disease_words))

    for p in papers:
        tl = (p.get("title") or "").lower()
        al = (p.get("abstract") or "").lower()

        off_t = official_l in tl
        off_a = official_l in al
        dis_t = _d_match(p.get("title") or "")
        dis_a = _d_match(p.get("abstract") or "")
        syn_t = any(s in tl for s in gene_syns_l)
        syn_a = any(s in al for s in gene_syns_l)

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
    """Europe PMC + PubMed E-utilities から論文を収集してスコア順に返す。"""

    # ── シノニム + MeSH ───────────────────────────────────────────────────────
    gene_syns    = _get_gene_synonyms(gene)
    mesh_heading = _get_mesh_heading(disease)
    # 疾患名の候補（元の名前 + MeSH 見出し語）
    disease_alts = list({disease, mesh_heading})

    print(f"    [文献] 遺伝子シノニム ({len(gene_syns)}件): "
          f"{', '.join(gene_syns[:4])}")
    print(f"    [文献] 疾患 MeSH: 「{mesh_heading}」")

    # ── Europe PMC 検索 ────────────────────────────────────────────────────────
    # EPMC はクォートなしの広いクエリで検索し EPMC 自身の関連度ランキングを活用。
    # プレプリント(PPR)・略称(CAH 等)を使う論文も取りこぼさないようにする。
    seen_ids: set[str] = set()
    epmc_papers: list[dict] = []

    def _add_epmc(items: list[dict]):
        for item in items:
            p = _epmc_to_paper(item)
            if p["pmid"] and p["pmid"] not in seen_ids:
                seen_ids.add(p["pmid"])
                epmc_papers.append(p)

    # Tier 1: クォートなしフリーテキスト（最広・プレプリント含む）
    _add_epmc(_epmc_search(f'{gene} {disease}', page_size=100, max_pages=4))

    # Tier 2: MeSH 別名でさらに補完
    if mesh_heading.lower() != disease.lower() and len(epmc_papers) < max_results * 2:
        _add_epmc(_epmc_search(f'{gene} {mesh_heading}', page_size=50, max_pages=2))

    # Tier 3: 遺伝子シノニム × 疾患名
    for syn in gene_syns[1:4]:
        if len(epmc_papers) >= max_results * 3:
            break
        _add_epmc(_epmc_search(f'{syn} {disease}', page_size=30, max_pages=1))

    print(f"    [文献] Europe PMC: {len(epmc_papers)} 件")

    # ── PubMed 補完（NCBI E-utilities） ─────────────────────────────────────
    supplement_pmids = _pubmed_supplement(
        gene, gene_syns, disease, mesh_heading,
        exclude_pmids=seen_ids,
        max_ids=max_results * 3,
    )
    time.sleep(0.3)
    pubmed_papers = _pubmed_fetch_papers(supplement_pmids[:300])
    print(f"    [文献] PubMed 補完: {len(pubmed_papers)} 件")

    # ── マージ + スコアリング ─────────────────────────────────────────────────
    all_papers = epmc_papers + pubmed_papers
    _score_papers(all_papers, gene, gene_syns, disease, disease_alts)

    # EPMC 論文: スコア関係なく保持（EPMC の関連度ランキングを信頼）
    # PubMed 補完: スコア 0（遺伝子名も疾患名も本文に出てこない）は除外
    pubmed_pmids = {p["pmid"] for p in pubmed_papers}
    all_papers = [
        p for p in all_papers
        if p["pmid"] not in pubmed_pmids or p["relevance_score"] > 0
    ]

    # ソート: 臨床研究 > スコア > 年 > EPMC優先
    all_papers.sort(
        key=lambda p: (
            p["is_clinical"],
            p["relevance_score"],
            p.get("year", ""),
            0 if p.get("_source") == "epmc" else 1,
        ),
        reverse=True,
    )

    # _source フィールドは内部用なので除去
    for p in all_papers:
        p.pop("_source", None)

    return all_papers[:max_results]
