"""GWAS Catalog REST API (EMBL-EBI, free for all use)."""
import requests

BASE = "https://www.ebi.ac.uk/gwas/rest/api"

def get_gwas_associations(gene_symbol: str, disease_query: str = None) -> list[dict]:
    """Return GWAS hits for a gene, optionally filtered by trait."""
    r = requests.get(f"{BASE}/genes/{gene_symbol}/associations", params={
        "projection": "associationByGene"
    }, timeout=20)

    if r.status_code == 404:
        return []
    r.raise_for_status()

    embedded = r.json().get("_embedded", {})
    associations = embedded.get("associations", [])

    results = []
    for assoc in associations[:20]:
        trait = assoc.get("efoTraits", [{}])[0].get("trait", "") if assoc.get("efoTraits") else ""
        pval = assoc.get("pvalue")
        or_beta = assoc.get("orPerCopyNum") or assoc.get("betaNum")
        snps = [s.get("rsId", "") for s in assoc.get("snps", [])]

        if disease_query and disease_query.lower() not in trait.lower():
            continue

        results.append({
            "trait": trait,
            "p_value": pval,
            "or_beta": or_beta,
            "snps": snps,
            "risk_allele_frequency": assoc.get("riskFrequency"),
        })

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
