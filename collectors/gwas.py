"""GWAS Catalog REST API (EMBL-EBI, free for all use)."""
import requests
from concurrent.futures import ThreadPoolExecutor

BASE = "https://www.ebi.ac.uk/gwas/rest/api"


def _fetch_snp_hits(snp: dict) -> list[dict]:
    """1 SNP の associations と trait 名を取得（並列実行用）。"""
    rsid = snp.get("rsId", "")
    assoc_link = (snp.get("_links", {}).get("associations") or {}).get("href")
    if not assoc_link:
        return []
    try:
        ra = requests.get(assoc_link, timeout=12)
        ra.raise_for_status()
        associations = ra.json().get("_embedded", {}).get("associations", [])
    except Exception:
        return []

    hits = []
    for assoc in associations:
        trait = ""
        tl = (assoc.get("_links", {}).get("efoTraits") or {}).get("href")
        if tl:
            try:
                rt = requests.get(tl, timeout=8)
                traits = rt.json().get("_embedded", {}).get("efoTraits", [])
                trait = ", ".join(t.get("trait", "") for t in traits if t.get("trait"))
            except Exception:
                pass
        hits.append({
            "trait":                 trait,
            "p_value":               assoc.get("pvalue"),
            "or_beta":               assoc.get("orPerCopyNum") or assoc.get("betaNum"),
            "snps":                  [rsid],
            "risk_allele_frequency": assoc.get("riskFrequency"),
        })
    return hits


def get_gwas_associations(gene_symbol: str, disease_query: str = None,
                          max_snps: int = 15, max_results: int = 20) -> list[dict]:
    """Return GWAS hits for a gene, optionally filtered by trait.

    旧 /genes/{gene}/associations は廃止 (500) されたため、
    findByGene で SNP を取得し、各 SNP の associations → efoTraits を並列で辿る。
    """
    try:
        r = requests.get(
            f"{BASE}/singleNucleotidePolymorphisms/search/findByGene",
            params={"geneName": gene_symbol, "size": max_snps}, timeout=20)
        if r.status_code == 404:
            return []
        r.raise_for_status()
        snps = r.json().get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
    except Exception:
        return []

    # SNP ごとの取得を並列化（逐次だと ~60s → 並列で数秒）
    # 語順違い（例: "type 2 diabetes mellitus" vs "diabetes mellitus type 2"）で
    # 単純な部分文字列一致だと本来ヒットすべき trait を取りこぼすため、
    # クエリを単語分割し全単語がトレイト文字列に含まれるかで判定する
    query_words = disease_query.lower().split() if disease_query else []

    results, seen = [], set()
    with ThreadPoolExecutor(max_workers=10) as ex:
        for hits in ex.map(_fetch_snp_hits, snps):
            for h in hits:
                trait_lower = (h["trait"] or "").lower()
                if query_words and not all(w in trait_lower for w in query_words):
                    continue
                key = (h["snps"][0], h["trait"])
                if key in seen:
                    continue
                seen.add(key)
                results.append(h)

    def _pv(x):
        try:
            return float(x.get("p_value") or 1)
        except (TypeError, ValueError):
            return 1.0
    results.sort(key=_pv)
    return results[:max_results]


def get_clinvar_variants(gene_symbol: str) -> list[dict]:
    """Return ClinVar pathogenic variants for a gene via NCBI API (public domain)."""
    import time
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    query = f"{gene_symbol}[Gene Name] AND (Pathogenic[Clinical significance] OR Likely pathogenic[Clinical significance])"

    r = requests.get(f"{base}/esearch.fcgi", params={
        "db": "clinvar", "term": query, "retmax": 100, "retmode": "json"
    }, timeout=15)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])

    if not ids:
        return []

    # 429対策: リトライ付きで esummary を呼ぶ
    for attempt in range(3):
        time.sleep(1 + attempt * 2)  # 1s, 3s, 5s
        r2 = requests.get(f"{base}/esummary.fcgi", params={
            "db": "clinvar", "id": ",".join(ids), "retmode": "json"
        }, timeout=15)
        if r2.status_code == 429:
            continue
        r2.raise_for_status()
        break
    else:
        return []  # リトライ上限に達した場合は空を返す
    result = r2.json().get("result", {})

    variants = []
    for vid in ids:
        if vid not in result:
            continue
        item = result[vid]
        # NCBI eutils は germline_classification（旧 clinical_significance）に
        # 臨床的意義・レビュー状態・trait を格納する
        cls = item.get("germline_classification") or item.get("clinical_significance") or {}
        traits = cls.get("trait_set") or []
        variants.append({
            "variant_id": vid,
            "title": item.get("title", ""),
            "clinical_significance": cls.get("description", ""),
            "condition": traits[0].get("trait_name", "") if traits else "",
            "review_status": cls.get("review_status", ""),
        })

    return variants
