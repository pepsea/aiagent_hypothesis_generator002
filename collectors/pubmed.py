"""PubMed literature collector via NCBI E-utilities (public domain).

シノニム対応の関連度スコアリング:
  スコア4: 公式シンボル AND 疾患名 → タイトルに一致
  スコア3: 公式シンボル AND 疾患名 → アブストラクトに一致
  スコア2: シノニム AND 疾患名 → タイトルに一致
  スコア1: シノニム AND 疾患名 → アブストラクトに一致
  スコア0: 遺伝子または疾患の片方のみ（内容参照）
新しい論文ほど優先（年順で2次ソート）。
"""
import time
import re
import requests
from typing import Optional

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get_gene_synonyms(gene_symbol: str) -> list[str]:
    """NCBI Gene API からシノニムを取得し、公式シンボル先頭のリストを返す。"""
    synonyms = [gene_symbol]

    try:
        # NCBI Gene 検索
        r = requests.get(f"{BASE}/esearch.fcgi", params={
            "db": "gene", "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
            "retmax": 1, "retmode": "json",
        }, timeout=10)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return synonyms

        time.sleep(0.3)

        # Gene エントリ取得
        r2 = requests.get(f"{BASE}/efetch.fcgi", params={
            "db": "gene", "id": ids[0], "rettype": "gene_table", "retmode": "text",
        }, timeout=10)
        r2.raise_for_status()

        # テキストからシノニム行を抽出
        for line in r2.text.splitlines():
            if line.startswith("Also known as"):
                raw = line.replace("Also known as", "").strip()
                for s in re.split(r"[;,]", raw):
                    s = s.strip()
                    if s and s.upper() != gene_symbol.upper():
                        synonyms.append(s)
                break

    except Exception:
        pass

    # UniProt からも補完（既存の uniprot モジュールを利用）
    try:
        from collectors.uniprot import get_protein_info
        info = get_protein_info(gene_symbol)
        # UniProt の keywords にシノニムが含まれる場合がある
        # gene_names フィールドを直接取りに行く
        r3 = requests.get("https://rest.uniprot.org/uniprotkb/search", params={
            "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
            "fields": "gene_names",
            "format": "json",
            "size": 1,
        }, timeout=10)
        r3.raise_for_status()
        results = r3.json().get("results", [])
        if results:
            genes = results[0].get("genes", [])
            for g in genes:
                # synonyms
                for syn in g.get("synonyms", []):
                    v = syn.get("value", "").strip()
                    if v and v not in synonyms:
                        synonyms.append(v)
                # ORF names
                for orf in g.get("orfNames", []):
                    v = orf.get("value", "").strip()
                    if v and v not in synonyms:
                        synonyms.append(v)
    except Exception:
        pass

    return synonyms


def search_pubmed(gene: str, disease: str, max_results: int = 10) -> list[dict]:
    """Return top PubMed abstracts scored by synonym-aware relevance."""

    # シノニム取得
    synonyms = get_gene_synonyms(gene)
    print(f"    遺伝子シノニム ({len(synonyms)}件): {', '.join(synonyms[:8])}{'...' if len(synonyms) > 8 else ''}")

    def _search(term: str, retmax: int) -> list[str]:
        r = requests.get(f"{BASE}/esearch.fcgi", params={
            "db": "pubmed", "term": term, "retmax": retmax,
            "retmode": "json", "sort": "pub+date",
        }, timeout=15)
        r.raise_for_status()
        return r.json().get("esearchresult", {}).get("idlist", [])

    # 公式シンボル検索
    queries_primary = [
        f'("{gene}"[Title/Abstract]) AND ("{disease}"[Title/Abstract])',
        f'"{gene}"[Gene Name] AND "{disease}"[MeSH Terms]',
        f"{gene} AND {disease}",
    ]
    # シノニム検索（上位3件まで）
    syn_terms = [f'"{s}"[Title/Abstract]' for s in synonyms[1:4] if len(s) >= 3]
    queries_synonym = []
    if syn_terms:
        syn_or = " OR ".join(syn_terms)
        queries_synonym = [
            f'({syn_or}) AND ("{disease}"[Title/Abstract])',
            f'({syn_or}) AND {disease}',
        ]

    # 収集（重複を除去しながら）
    seen = set()
    all_ids = []

    for q in queries_primary + queries_synonym:
        ids = _search(q, retmax=max_results * 2)
        for i in ids:
            if i not in seen:
                seen.add(i)
                all_ids.append(i)
        if len(all_ids) >= max_results * 3:
            break

    if not all_ids:
        return []

    time.sleep(0.5)

    # サマリー取得
    r = requests.post(f"{BASE}/esummary.fcgi", data={
        "db": "pubmed", "id": ",".join(all_ids[:40]), "retmode": "json"
    }, timeout=15)
    r.raise_for_status()
    result = r.json().get("result", {})

    papers = []
    for pmid in all_ids[:40]:
        if pmid not in result:
            continue
        item = result[pmid]
        papers.append({
            "pmid": pmid,
            "title": item.get("title", ""),
            "journal": item.get("fulljournalname", ""),
            "year": item.get("pubdate", "")[:4],
            "authors": [a.get("name", "") for a in item.get("authors", [])[:3]],
            "abstract": "",
            "relevance_score": 0,
            "match_type": "",
        })

    # アブストラクト取得 & スコアリング
    _score_papers(papers, gene, synonyms, disease)

    # score 降順 → year 降順
    papers.sort(key=lambda p: (p["relevance_score"], p["year"]), reverse=True)

    return papers[:max_results]


def _score_papers(papers: list[dict], gene: str, synonyms: list[str], disease: str):
    """Fetch abstracts and compute synonym-aware relevance scores in-place."""
    if not papers:
        return

    pmids = [p["pmid"] for p in papers]
    time.sleep(0.3)

    try:
        r = requests.get(f"{BASE}/efetch.fcgi", params={
            "db": "pubmed", "id": ",".join(pmids),
            "rettype": "abstract", "retmode": "xml",
        }, timeout=25)
        r.raise_for_status()
        xml = r.text
    except Exception:
        return

    # 簡易XMLパース
    abstract_blocks = re.findall(
        r"<MedlineCitation[^>]*>.*?<PMID[^>]*>(\d+)</PMID>.*?<AbstractText[^>]*>(.*?)</AbstractText>.*?</MedlineCitation>",
        xml, re.DOTALL
    )
    abstract_map = {pmid: re.sub(r"<[^>]+>", " ", text).strip()
                    for pmid, text in abstract_blocks}

    official_l = gene.lower()
    synonyms_l = [s.lower() for s in synonyms[1:]]  # 公式シンボルを除くシノニム
    disease_words = [w for w in disease.lower().split() if len(w) > 4]

    def disease_match(text: str) -> bool:
        t = text.lower()
        return disease.lower() in t or all(w in t for w in disease_words)

    for p in papers:
        abstract = abstract_map.get(p["pmid"], "")
        p["abstract"] = abstract

        title_l    = p["title"].lower()
        abstract_l = abstract.lower()

        official_in_title = official_l in title_l
        official_in_abs   = official_l in abstract_l
        disease_in_title  = disease_match(p["title"])
        disease_in_abs    = disease_match(abstract)
        syn_in_title      = any(s in title_l for s in synonyms_l)
        syn_in_abs        = any(s in abstract_l for s in synonyms_l)

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
        p["match_type"] = match


def fetch_abstract(pmid: str) -> Optional[str]:
    r = requests.get(f"{BASE}/efetch.fcgi", params={
        "db": "pubmed", "id": pmid, "rettype": "abstract", "retmode": "text"
    }, timeout=15)
    r.raise_for_status()
    return r.text.strip()
