"""HPO (Human Phenotype Ontology) API collector.

Retrieves disease phenotypes and their associated genes via hpo.jax.org API,
then evaluates overlap between a target gene's PPI partners and HPO-gene sets.
"""
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

HPO_BASE = "https://hpo.jax.org/api/hpo"
_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json"})


def _get(url: str, params: dict = None, timeout: int = 15) -> dict:
    r = _SESSION.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.json()


def search_disease(disease_name: str, max_results: int = 5) -> list[dict]:
    """Search for diseases by name. Returns list of {diseaseId, diseaseName}."""
    data = _get(f"{HPO_BASE}/search", params={
        "q": disease_name, "max": max_results, "category": "diseases"
    })
    return data.get("diseases") or []


def get_disease_phenotypes(disease_id: str) -> dict:
    """Get phenotype terms and associated genes for a disease.

    disease_id: e.g. "OMIM:201910" or "ORPHA:12345"
    Returns {
        "disease_id": str,
        "disease_name": str,
        "phenotypes": [{"hpo_id", "name", "frequency"}, ...],
        "genes": [{"gene_symbol", "gene_id"}, ...],
    }
    """
    data = _get(f"{HPO_BASE}/disease/{disease_id}")
    phenotypes_raw = data.get("catTermsCombo") or data.get("phenotypes") or []
    phenotypes = []
    for p in phenotypes_raw:
        hpo_id = p.get("ontologyId") or (p.get("term") or {}).get("id")
        name   = p.get("name") or (p.get("term") or {}).get("name")
        freq   = p.get("frequency") or {}
        freq_label = freq.get("label") if isinstance(freq, dict) else str(freq or "")
        if hpo_id and name:
            phenotypes.append({"hpo_id": hpo_id, "name": name, "frequency": freq_label})

    genes_raw = data.get("associatedGenes") or data.get("geneAssoc") or []
    genes = []
    for g in genes_raw:
        sym = g.get("geneSymbol") or g.get("symbol") or g.get("gene")
        gid = g.get("geneId") or g.get("entrezId")
        if sym:
            genes.append({"gene_symbol": sym, "gene_id": gid})

    return {
        "disease_id":   disease_id,
        "disease_name": data.get("disease", {}).get("diseaseName", "") if isinstance(data.get("disease"), dict) else "",
        "phenotypes":   phenotypes,
        "genes":        genes,
    }


def get_term_genes(hpo_id: str, max_results: int = 200) -> list[str]:
    """Get gene symbols associated with a specific HPO term."""
    try:
        data = _get(f"{HPO_BASE}/term/{hpo_id}/genes",
                    params={"max": max_results})
        genes = data.get("genes") or data.get("geneAssoc") or []
        return [
            g.get("geneSymbol") or g.get("symbol") or g.get("gene", "")
            for g in genes
            if g.get("geneSymbol") or g.get("symbol") or g.get("gene")
        ]
    except Exception:
        return []


def evaluate_ppi_hpo_overlap(
    gene: str,
    ppi_partners: list[str],
    disease_name: str,
    omim_id: str = None,
    max_phenotypes: int = 30,
) -> dict:
    """Evaluate overlap between PPI partners and HPO symptom-associated genes.

    Args:
        gene: target gene symbol
        ppi_partners: list of PPI partner gene symbols (hub-filtered)
        disease_name: disease name for HPO lookup
        omim_id: optional OMIM ID (e.g. "OMIM:201910") to skip disease search
        max_phenotypes: maximum number of HPO terms to query

    Returns dict with:
        disease_id, disease_name, hpo_term_count,
        per_term: [{hpo_id, name, frequency, overlap_genes, overlap_count}],
        summary: {total_hpo_genes, overlap_genes, overlap_count, overlap_score,
                  top_genes [{symbol, term_count}]}
    """
    ppi_set = {p.upper() for p in ppi_partners if p}
    gene_upper = gene.upper()

    # Step 1: resolve disease to HPO disease ID
    disease_id = omim_id
    disease_label = disease_name
    if not disease_id:
        hits = search_disease(disease_name, max_results=3)
        if hits:
            disease_id    = hits[0].get("diseaseId")
            disease_label = hits[0].get("diseaseName", disease_name)

    if not disease_id:
        return {"error": "HPO disease not found", "disease_name": disease_name}

    # Step 2: get phenotypes for the disease
    try:
        disease_data = get_disease_phenotypes(disease_id)
    except Exception as e:
        return {"error": f"HPO API error: {e}", "disease_id": disease_id}

    phenotypes = disease_data.get("phenotypes", [])[:max_phenotypes]
    disease_label = disease_data.get("disease_name") or disease_label

    # Step 3: for each phenotype term, get associated genes in parallel
    def _fetch_term(pheno):
        hpo_id = pheno["hpo_id"]
        genes  = get_term_genes(hpo_id, max_results=300)
        return hpo_id, genes

    term_genes: dict[str, list[str]] = {}
    if phenotypes:
        with ThreadPoolExecutor(max_workers=8) as ex:
            futs = {ex.submit(_fetch_term, p): p for p in phenotypes}
            for fut in as_completed(futs):
                try:
                    hpo_id, genes = fut.result()
                    term_genes[hpo_id] = genes
                except Exception:
                    pass

    # Step 4: calculate per-term overlap with PPI partners
    per_term = []
    all_hpo_genes: set[str] = set()
    for pheno in phenotypes:
        hpo_id = pheno["hpo_id"]
        tgenes = term_genes.get(hpo_id, [])
        tgenes_upper = {g.upper() for g in tgenes if g}
        all_hpo_genes.update(tgenes_upper)
        overlap = sorted(ppi_set & tgenes_upper)
        per_term.append({
            "hpo_id":       hpo_id,
            "name":         pheno["name"],
            "frequency":    pheno.get("frequency", ""),
            "hpo_gene_count": len(tgenes_upper),
            "overlap_genes": [g for g in overlap],
            "overlap_count": len(overlap),
        })

    per_term.sort(key=lambda x: x["overlap_count"], reverse=True)

    # Step 5: aggregate statistics
    overlap_genes = ppi_set & all_hpo_genes
    # gene frequency across terms
    gene_term_count: dict[str, int] = {}
    for pt in per_term:
        for g in pt["overlap_genes"]:
            gene_term_count[g.upper()] = gene_term_count.get(g.upper(), 0) + 1

    top_genes = sorted(
        [{"symbol": g, "term_count": c} for g, c in gene_term_count.items()],
        key=lambda x: x["term_count"], reverse=True
    )[:15]

    total_hpo_genes = len(all_hpo_genes)
    overlap_score = round(
        len(overlap_genes) / max(1, min(len(ppi_set), total_hpo_genes)), 3
    ) if (ppi_set and total_hpo_genes) else 0.0

    # Also check if target gene itself is in HPO gene set
    target_in_hpo = gene_upper in all_hpo_genes

    return {
        "disease_id":       disease_id,
        "disease_name":     disease_label,
        "hpo_term_count":   len(phenotypes),
        "per_term":         per_term,
        "summary": {
            "total_hpo_genes":  total_hpo_genes,
            "ppi_partner_count": len(ppi_set),
            "overlap_genes":    sorted(list(overlap_genes)),
            "overlap_count":    len(overlap_genes),
            "overlap_score":    overlap_score,
            "target_in_hpo":    target_in_hpo,
            "top_genes":        top_genes,
        },
    }
