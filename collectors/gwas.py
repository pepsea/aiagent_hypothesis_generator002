"""GWAS Catalog REST API (EMBL-EBI, free for all use)."""
import requests

BASE = "https://www.ebi.ac.uk/gwas/rest/api"

def get_gwas_associations(gene_symbol: str, disease_query: str = None,
                          max_snps: int = 25, max_results: int = 20) -> list[dict]:
    """Return GWAS hits for a gene, optionally filtered by trait.

    旧 /genes/{gene}/associations は廃止 (500) されたため、
    findByGene で SNP を取得し、各 SNP の associations → efoTraits を辿る。
    """
    # 1. 遺伝子にマップされる SNP を取得
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

    results = []
    seen = set()
    for snp in snps:
        if len(results) >= max_results:
            break
        rsid = snp.get("rsId", "")
        assoc_link = (snp.get("_links", {}).get("associations") or {}).get("href")
        if not assoc_link:
            continue
        try:
            ra = requests.get(assoc_link, timeout=15)
            ra.raise_for_status()
            associations = ra.json().get("_embedded", {}).get("associations", [])
        except Exception:
            continue

        for assoc in associations:
            pval = assoc.get("pvalue")
            or_beta = assoc.get("orPerCopyNum") or assoc.get("betaNum")

            # trait 名を efoTraits リンクから取得
            trait = ""
            tl = (assoc.get("_links", {}).get("efoTraits") or {}).get("href")
            if tl:
                try:
                    rt = requests.get(tl, timeout=10)
                    traits = rt.json().get("_embedded", {}).get("efoTraits", [])
                    trait = ", ".join(t.get("trait", "") for t in traits if t.get("trait"))
                except Exception:
                    pass

            if disease_query and disease_query.lower() not in trait.lower():
                continue

            key = (rsid, trait)
            if key in seen:
                continue
            seen.add(key)

            results.append({
                "trait":                 trait,
                "p_value":               pval,
                "or_beta":               or_beta,
                "snps":                  [rsid],
                "risk_allele_frequency": assoc.get("riskFrequency"),
            })
            if len(results) >= max_results:
                break

    # p値昇順（数値化できるもの優先）
    def _pv(x):
        try:
            return float(x.get("p_value") or 1)
        except (TypeError, ValueError):
            return 1.0
    results.sort(key=_pv)
    return results


def get_clinvar_variants(gene_symbol: str) -> list[dict]:
    """Return ClinVar pathogenic variants for a gene via NCBI API (public domain)."""
    import time
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    query = f"{gene_symbol}[Gene Name] AND (Pathogenic[Clinical significance] OR Likely pathogenic[Clinical significance])"

    r = requests.get(f"{base}/esearch.fcgi", params={
        "db": "clinvar", "term": query, "retmax": 10, "retmode": "json"
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
        variants.append({
            "variant_id": vid,
            "title": item.get("title", ""),
            "clinical_significance": item.get("clinical_significance", {}).get("description", ""),
            "condition": item.get("trait_set", [{}])[0].get("trait_name", "") if item.get("trait_set") else "",
            "review_status": item.get("clinical_significance", {}).get("review_status", ""),
        })

    return variants
