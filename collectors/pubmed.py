"""PubMed literature collector via NCBI E-utilities (public domain).

検索戦略（論文漏れを最小化）:
  Tier 1: gene[Gene/Protein Name] AND disease[MeSH Terms]       ← MeSH 自動マッピングを活用
  Tier 2: "gene"[Title/Abstract] AND disease[MeSH Terms]        ← MeSH インデックス済み論文
  Tier 3: "gene"[Title/Abstract] AND "disease"[Title/Abstract]  ← 未インデックス/新着論文
  Tier 4: gene_synonyms[Title/Abstract] AND disease[MeSH Terms] ← 旧称/別称での論文
  Tier 5: gene_synonyms AND disease_synonyms[Title/Abstract]    ← 略称・別名の組み合わせ

スコアリング（関連度）:
  4: 公式シンボル × 疾患名 → タイトル一致
  3: 公式シンボル × 疾患名 → アブストラクト一致
  2: シノニム × 疾患名 → タイトル一致
  1: シノニム × 疾患名 → アブストラクト一致
  0: 片方のみ（内容参照）

新しい論文ほど優先（年順で2次ソート）。
"""
import time
import re
import requests
from typing import Optional

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# PubMed Publication Type のうち「臨床（ヒトでの研究）」とみなすもの。
# これに該当しない論文（in vitro/動物モデル/総説/基礎研究等）は非臨床として扱う。
# "Comparative Study" は動物実験の群間比較にも routine に付与されるため対象外
# （ヒト臨床研究であることを保証しない）。"Meta-Analysis" も対象論文の種類を
# 保証しないため除外。
CLINICAL_PUBTYPES = {
    "clinical trial", "clinical trial, phase i", "clinical trial, phase ii",
    "clinical trial, phase iii", "clinical trial, phase iv",
    "controlled clinical trial", "randomized controlled trial",
    "pragmatic clinical trial", "adaptive clinical trial",
    "observational study", "multicenter study", "case reports",
}


def _is_clinical(pub_types: list[str]) -> bool:
    """pubtype リストに臨床研究を示す種別が含まれるか判定する。"""
    return any((pt or "").lower() in CLINICAL_PUBTYPES for pt in (pub_types or []))

# MeSH クエリ用ヘッダー（MeSH の見出し語順 "Dystrophy, Duchenne Muscular" など）を保持
_mesh_heading_cache: dict[str, str] = {}


# ─────────────────────────────────────────────────────────────────────────────
# 遺伝子シノニム取得
# ─────────────────────────────────────────────────────────────────────────────

def get_gene_synonyms(gene_symbol: str) -> list[str]:
    """NCBI Gene + UniProt から遺伝子シノニムを取得。公式シンボルを先頭に返す。"""
    synonyms = [gene_symbol]
    seen = {gene_symbol.upper()}

    def add(name: str):
        name = name.strip()
        if name and name.upper() not in seen and len(name) >= 2:
            seen.add(name.upper())
            synonyms.append(name)

    # ── NCBI Gene ──────────────────────────────────────────────────────────
    try:
        r = requests.get(f"{BASE}/esearch.fcgi", params={
            "db": "gene",
            "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
            "retmax": 1, "retmode": "json",
        }, timeout=10)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if ids:
            time.sleep(0.3)
            r2 = requests.get(f"{BASE}/efetch.fcgi", params={
                "db": "gene", "id": ids[0],
                "rettype": "gene_table", "retmode": "text",
            }, timeout=10)
            r2.raise_for_status()
            for line in r2.text.splitlines():
                if line.startswith("Also known as"):
                    raw = line.replace("Also known as", "").strip()
                    for s in re.split(r"[;,]", raw):
                        add(s)
                    break
    except Exception:
        pass

    # ── UniProt gene_names + protein_name ──────────────────────────────────
    try:
        r3 = requests.get("https://rest.uniprot.org/uniprotkb/search", params={
            "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
            "fields": "gene_names,protein_name",
            "format": "json",
            "size": 1,
        }, timeout=10)
        r3.raise_for_status()
        results = r3.json().get("results", [])
        if results:
            entry = results[0]
            # 遺伝子シノニム
            for g in entry.get("genes", []):
                for syn in g.get("synonyms", []):
                    add(syn.get("value", ""))
                for orf in g.get("orfNames", []):
                    add(orf.get("value", ""))
            # タンパク質名（推奨名・別名・略称）
            pd = entry.get("proteinDescription", {})
            rec = pd.get("recommendedName", {})
            add(rec.get("fullName", {}).get("value", ""))
            for sn in rec.get("shortNames", []):
                add(sn.get("value", ""))
            for alt in pd.get("alternativeNames", []):
                add(alt.get("fullName", {}).get("value", ""))
                for sn in alt.get("shortNames", []):
                    add(sn.get("value", ""))
    except Exception:
        pass

    return synonyms


# ─────────────────────────────────────────────────────────────────────────────
# 疾患シノニム取得
# ─────────────────────────────────────────────────────────────────────────────

def get_disease_synonyms(disease: str, efo_id: str = None) -> tuple[list[str], str]:
    """疾患名のシノニムと MeSH heading を返す。

    Returns
    -------
    synonyms : list[str]   公式名を先頭に含む別名リスト（略称含む）
    mesh_heading : str     PubMed [MeSH Terms] クエリに使う見出し語（空文字の場合あり）
    """
    synonyms = [disease]
    seen = {disease.lower()}
    mesh_heading = ""

    def add(name: str):
        name = name.strip()
        if name and name.lower() not in seen and len(name) >= 2:
            seen.add(name.lower())
            synonyms.append(name)

    # ── 1. NCBI MeSH entry terms ────────────────────────────────────────────
    # efetch?db=mesh はプレーンテキスト形式で返る
    try:
        # "[MeSH Subheading]" クエリが MeSH descriptor に正確にマップされる
        for mesh_query in [f"{disease}[MeSH Subheading]", disease]:
            r = requests.get(f"{BASE}/esearch.fcgi", params={
                "db": "mesh", "term": mesh_query,
                "retmax": 1, "retmode": "json",
            }, timeout=10)
            r.raise_for_status()
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                break

        if ids:
            time.sleep(0.3)
            r2 = requests.get(f"{BASE}/efetch.fcgi", params={
                "db": "mesh", "id": ids[0],
            }, timeout=10)
            r2.raise_for_status()
            text = r2.text

            lines = text.splitlines()

            # "N: <DescriptorName>" 行を探す → MeSH heading
            for line in lines:
                m = re.match(r"^\d+:\s+(.+)$", line.strip())
                if m:
                    mesh_heading = m.group(1).strip()
                    if mesh_heading.lower() != disease.lower():
                        add(mesh_heading)
                    break

            # "Entry Terms:" セクションを解析
            # ツリー階層（"All MeSH Categories" 以降）は空行で区切られているので
            # 最初の空行で終了する
            in_entry = False
            for line in lines:
                if line.strip().startswith("Entry Terms:"):
                    in_entry = True
                    continue
                if in_entry:
                    if line.strip() == "":
                        break   # 空行でentry terms終了（ツリー階層の手前）
                    term = line.strip()
                    if term:
                        add(term)
    except Exception:
        pass

    # ── 2. OLS4 / EFO・MONDO synonyms ──────────────────────────────────────
    if efo_id:
        try:
            if efo_id.startswith("EFO_"):
                iri  = f"http://www.ebi.ac.uk/efo/{efo_id}"
                onto = "efo"
            elif efo_id.startswith("MONDO_"):
                iri  = f"http://purl.obolibrary.org/obo/{efo_id.replace(':', '_')}"
                onto = "mondo"
            elif efo_id.startswith("Orphanet_") or efo_id.startswith("HP_"):
                iri  = f"http://www.orpha.net/ORDO/{efo_id}"
                onto = "ordo"
            else:
                iri = None

            if iri:
                r3 = requests.get(
                    f"https://www.ebi.ac.uk/ols4/api/ontologies/{onto}/terms",
                    params={"iri": iri}, timeout=12,
                )
                r3.raise_for_status()
                terms_list = (r3.json().get("_embedded") or {}).get("terms", [])
                if terms_list:
                    t = terms_list[0]
                    add(t.get("label", ""))
                    for s in t.get("synonyms") or []:
                        add(s)
                    for oa in t.get("obo_synonym") or []:
                        add(oa.get("name", ""))
        except Exception:
            pass

    return synonyms, mesh_heading


# ─────────────────────────────────────────────────────────────────────────────
# PubMed 検索
# ─────────────────────────────────────────────────────────────────────────────

def search_pubmed(
    gene: str,
    disease: str,
    max_results: int = 100,
    disease_efo_id: str = None,
) -> list[dict]:
    """Return top PubMed abstracts scored by synonym-aware relevance.

    Parameters
    ----------
    gene            : HGNC gene symbol
    disease         : disease name (from OpenTargets or user input)
    max_results     : number of papers to return (after scoring). Default
                      100 — display/report should note this is a "top 100,
                      most recent" list, not an exhaustive one.
    disease_efo_id  : EFO/MONDO ID for richer disease synonym lookup

    論文抽出は「疾患名と遺伝子名（シノニム含む）がタイトルまたはアブストラクトに
    記載されている」ものを中心とする。score==0（MeSH等の間接一致のみで、
    タイトル/アブストラクトに実際の語が出てこない「内容参照」papers）は除外する。
    """

    # ── シノニム取得 ─────────────────────────────────────────────────────────
    gene_syns = get_gene_synonyms(gene)
    disease_syns, mesh_heading = get_disease_synonyms(disease, efo_id=disease_efo_id)

    n_gsyn = len(gene_syns) - 1
    n_dsyn = len(disease_syns) - 1
    print(f"    遺伝子シノニム ({len(gene_syns)}件): "
          f"{', '.join(gene_syns[:6])}{'...' if n_gsyn > 5 else ''}")
    print(f"    疾患シノニム  ({len(disease_syns)}件): "
          f"{', '.join(disease_syns[:6])}{'...' if n_dsyn > 5 else ''}"
          + (f"  MeSH: 「{mesh_heading}」" if mesh_heading else ""))

    # ── クエリ構築ヘルパー ────────────────────────────────────────────────────
    def _q_term(term: str, field: str = "Title/Abstract") -> str:
        """1単語は完全一致。複数単語は各単語を AND 結合（部分一致）。
        3文字以下の短い語（略称など）は除外しない（PubMed が自動処理する）。
        """
        words = term.split()
        if len(words) <= 1:
            return f'"{term}"[{field}]'
        return "(" + " AND ".join(f'"{w}"[{field}]' for w in words) + ")"

    def _q_gene_official() -> str:
        return _q_term(gene)

    def _q_gene_syns(max_syn: int = 8) -> str:
        """公式シンボル + シノニムの OR クエリ（Title/Abstract）。"""
        terms = [_q_term(s) for s in gene_syns[:max_syn + 1]]
        return "(" + " OR ".join(terms) + ")"

    def _q_disease_mesh() -> str:
        heading = mesh_heading or disease
        return f'"{heading}"[MeSH Terms]'

    def _q_disease_text(max_syn: int = 4) -> str:
        """疾患名 + 略称シノニムの OR クエリ（Title/Abstract）。"""
        terms = [_q_term(s) for s in disease_syns[:max_syn + 1]]
        return "(" + " OR ".join(terms) + ")"

    def _q_clinical_types() -> str:
        return "(" + " OR ".join(f'"{t}"[Publication Type]' for t in [
            "Clinical Trial", "Clinical Trial, Phase I", "Clinical Trial, Phase II",
            "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
            "Controlled Clinical Trial", "Randomized Controlled Trial",
            "Observational Study", "Multicenter Study", "Case Reports",
        ]) + ")"

    def _search(term: str) -> list[str]:
        try:
            r = requests.get(f"{BASE}/esearch.fcgi", params={
                "db": "pubmed", "term": term,
                "retmax": max_results * 3,
                "retmode": "json", "sort": "pub+date",
            }, timeout=15)
            r.raise_for_status()
            return r.json().get("esearchresult", {}).get("idlist", [])
        except Exception:
            return []

    # ── 段階的クエリ（Tier 0 → 5） ───────────────────────────────────────────
    queries = [
        # Tier 0: 臨床研究を明示的に検索する。
        # 他のTierは全て sort=pub+date で「直近の論文」を優先的に取得するため、
        # 対象が昔に活発だった薬剤（例: BACE1阻害薬治験は2016-2019年が中心）の場合、
        # 近年の基礎研究論文に押し出され、臨床試験論文がそもそも候補プールに
        # 入らないことがある。臨床研究に限定した検索を別立てで行うことで、
        # 発表時期に関わらず臨床論文を確実に候補へ含める。
        f'{_q_gene_syns(max_syn=8)} AND {_q_disease_text(max_syn=4)} AND {_q_clinical_types()}',
        # Tier 1: 遺伝子フィールド × MeSH（最も精度が高い）
        f'"{gene}"[Gene/Protein Name] AND {_q_disease_mesh()}',
        # Tier 2: 公式シンボル × MeSH（新着論文もカバー）
        f'{_q_gene_official()} AND {_q_disease_mesh()}',
        # Tier 3: 公式シンボル × 疾患名テキスト（MeSH 未インデックス・新着）
        f'{_q_gene_official()} AND {_q_disease_text(max_syn=3)}',
        # Tier 4: 遺伝子シノニム × MeSH
        f'{_q_gene_syns(max_syn=8)} AND {_q_disease_mesh()}',
        # Tier 5: 遺伝子シノニム × 疾患テキスト（シノニム×略称の組み合わせ）
        f'{_q_gene_syns(max_syn=8)} AND {_q_disease_text(max_syn=4)}',
    ]

    seen: set[str] = set()
    all_ids: list[str] = []

    for q in queries:
        ids = _search(q)
        for i in ids:
            if i not in seen:
                seen.add(i)
                all_ids.append(i)
        if len(all_ids) >= max_results * 5:
            break
        time.sleep(0.2)

    if not all_ids:
        return []

    time.sleep(0.4)

    # ── サマリー取得 ──────────────────────────────────────────────────────────
    fetch_ids = all_ids[:min(max_results * 3, 300)]
    r = requests.post(f"{BASE}/esummary.fcgi", data={
        "db": "pubmed", "id": ",".join(fetch_ids), "retmode": "json",
    }, timeout=20)
    r.raise_for_status()
    result = r.json().get("result", {})

    papers = []
    for pmid in fetch_ids:
        if pmid not in result:
            continue
        item = result[pmid]
        pub_types = item.get("pubtype") or []
        papers.append({
            "pmid":            pmid,
            "title":           item.get("title", ""),
            "journal":         item.get("fulljournalname", ""),
            "year":            item.get("pubdate", "")[:4],
            "authors":         [a.get("name", "") for a in item.get("authors", [])[:3]],
            "abstract":        "",
            "relevance_score": 0,
            "match_type":      "",
            "pub_types":       pub_types,
            "is_clinical":     _is_clinical(pub_types),
        })

    # ── アブストラクト取得 & スコアリング ─────────────────────────────────────
    _score_papers(papers, gene, gene_syns, disease, disease_syns)

    # score==0（タイトル/アブストラクトに遺伝子名・疾患名いずれかの実際の語が
    # 出てこない「内容参照」papers、MeSH等の間接一致のみ）は常に除外する。
    # 論文が少ないペアの補強はパスウェイ隣接遺伝子の related_gene_papers で行う。
    papers = [p for p in papers if p["relevance_score"] > 0]

    # 臨床研究を優先し（is_clinical 降順）、その中で関連度→年の順に並べる。
    papers.sort(key=lambda p: (p["is_clinical"], p["relevance_score"], p["year"]), reverse=True)

    return papers[:max_results]


# ─────────────────────────────────────────────────────────────────────────────
# スコアリング
# ─────────────────────────────────────────────────────────────────────────────

def _score_papers(
    papers: list[dict],
    gene: str,
    gene_syns: list[str],
    disease: str,
    disease_syns: list[str],
):
    """Fetch abstracts and compute synonym-aware relevance scores in-place."""
    if not papers:
        return

    pmids = [p["pmid"] for p in papers]
    time.sleep(0.3)

    try:
        r = requests.get(f"{BASE}/efetch.fcgi", params={
            "db": "pubmed", "id": ",".join(pmids),
            "rettype": "abstract", "retmode": "xml",
        }, timeout=30)
        r.raise_for_status()
        xml = r.text
    except Exception:
        return

    # 簡易 XML パース
    abstract_map = {
        pmid: re.sub(r"<[^>]+>", " ", text).strip()
        for pmid, text in re.findall(
            r"<MedlineCitation[^>]*>.*?<PMID[^>]*>(\d+)</PMID>"
            r".*?<AbstractText[^>]*>(.*?)</AbstractText>"
            r".*?</MedlineCitation>",
            xml, re.DOTALL,
        )
    }

    official_l   = gene.lower()
    gene_syns_l  = [s.lower() for s in gene_syns[1:]]   # 公式シンボル除くシノニム
    disease_l    = disease.lower()
    # 疾患マッチ用: 全シノニムを使う（短い略称は誤ヒットしやすいが後段スコアで補正）
    disease_all_l = [s.lower() for s in disease_syns]
    disease_words = [w for w in disease_l.split() if len(w) > 4]

    def disease_match(text: str) -> bool:
        t = text.lower()
        return any(d in t for d in disease_all_l) or all(w in t for w in disease_words)

    for p in papers:
        abstract = abstract_map.get(p["pmid"], "")
        p["abstract"] = abstract

        title_l    = p["title"].lower()
        abstract_l = abstract.lower()

        official_in_title = official_l in title_l
        official_in_abs   = official_l in abstract_l
        disease_in_title  = disease_match(p["title"])
        disease_in_abs    = disease_match(abstract)
        syn_in_title      = any(s in title_l for s in gene_syns_l)
        syn_in_abs        = any(s in abstract_l for s in gene_syns_l)

        if official_in_title and disease_in_title:
            score, match = 4, "公式シンボル×タイトル"
        elif (official_in_title or official_in_abs) and (disease_in_title or disease_in_abs):
            score, match = 3, "公式シンボル×アブスト"
        elif syn_in_title and disease_in_title:
            score, match = 2, "シノニム×タイトル"
        elif (syn_in_title or syn_in_abs) and (disease_in_title or disease_in_abs):
            score, match = 1, "シノニム×アブスト"
        else:
            score, match = 0, "内容参照"

        p["relevance_score"] = score
        p["match_type"]      = match


# ─────────────────────────────────────────────────────────────────────────────
# ユーティリティ
# ─────────────────────────────────────────────────────────────────────────────

def fetch_abstract(pmid: str) -> Optional[str]:
    r = requests.get(f"{BASE}/efetch.fcgi", params={
        "db": "pubmed", "id": pmid,
        "rettype": "abstract", "retmode": "text",
    }, timeout=15)
    r.raise_for_status()
    return r.text.strip()
