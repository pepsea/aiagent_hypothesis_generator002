"""HPO (Human Phenotype Ontology) collector.

Disease phenotypes and associated genes via multiple sources:
  1. HPO JAX API  (hpo.jax.org/api/hpo/disease/{id})  — OMIM/ORPHA IDs
  2. Monarch Initiative API  (api.monarchinitiative.org) — MONDO IDs

Then evaluates overlap between a target gene's PPI partners and HPO-gene sets.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

HPO_BASE     = "https://hpo.jax.org/api/hpo"
MONARCH_BASE = "https://api.monarchinitiative.org/v3/api"

_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json"})


def _get(url: str, params: dict = None, timeout: int = 15) -> dict:
    r = _SESSION.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.json()


# ── Disease ID helpers ─────────────────────────────────────────────────────

def _mondo_to_hpo_format(mondo_id: str) -> str | None:
    """Convert MONDO_0008728 or MONDO:0008728 → OMIM ID via Monarch API.

    Returns "OMIM:201910" style string, or None if not found.
    """
    mid = mondo_id.replace("_", ":") if "_" in mondo_id else mondo_id
    if not mid.upper().startswith("MONDO"):
        return None
    try:
        data = _get(f"{MONARCH_BASE}/entity/{mid}")
        for xref in (data.get("xrefs") or []):
            if isinstance(xref, str) and xref.startswith("OMIM:"):
                return xref
            if isinstance(xref, dict):
                val = xref.get("id") or xref.get("curie", "")
                if val.startswith("OMIM:"):
                    return val
    except Exception:
        pass
    return None


# ── Disease → phenotypes ───────────────────────────────────────────────────

def _hpo_api_disease_phenotypes(disease_id: str) -> dict | None:
    """Try HPO JAX API for OMIM/ORPHA disease IDs."""
    try:
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
        genes = [
            {"gene_symbol": g.get("geneSymbol") or g.get("symbol") or g.get("gene"),
             "gene_id": g.get("geneId") or g.get("entrezId")}
            for g in genes_raw
            if g.get("geneSymbol") or g.get("symbol") or g.get("gene")
        ]
        dis = data.get("disease") or {}
        return {
            "disease_id":   disease_id,
            "disease_name": dis.get("diseaseName", "") if isinstance(dis, dict) else "",
            "phenotypes":   phenotypes,
            "genes":        genes,
        }
    except Exception:
        return None


def _monarch_disease_phenotypes(mondo_id: str) -> dict | None:
    """Use Monarch Initiative API for MONDO disease IDs."""
    mid = mondo_id.replace("_", ":") if "_" in mondo_id else mondo_id
    try:
        data = _get(f"{MONARCH_BASE}/association/disease/{mid}/phenotype",
                    params={"limit": 200})
        items = data.get("items") or data.get("associations") or []
        phenotypes = []
        for item in items:
            subj = item.get("object") or item.get("phenotype") or {}
            if isinstance(subj, str):
                hpo_id = subj
                name = subj
            else:
                hpo_id = subj.get("id") or subj.get("curie") or subj.get("identifier", "")
                name   = subj.get("label") or subj.get("name", "")
            freq = (item.get("frequency") or {})
            freq_label = freq.get("label", "") if isinstance(freq, dict) else ""
            if hpo_id and name:
                phenotypes.append({"hpo_id": hpo_id, "name": name, "frequency": freq_label})
        return {
            "disease_id":   mid,
            "disease_name": "",
            "phenotypes":   phenotypes,
            "genes":        [],
        }
    except Exception:
        return None


def get_disease_phenotypes(disease_id: str) -> dict:
    """Get HPO phenotypes for a disease. Accepts OMIM:, ORPHA:, MONDO: or MONDO_ IDs."""
    uid = disease_id.replace("_", ":") if "_" in disease_id else disease_id

    # MONDO ID → try to find OMIM equivalent first, then Monarch
    if uid.upper().startswith("MONDO"):
        omim_id = _mondo_to_hpo_format(uid)
        if omim_id:
            result = _hpo_api_disease_phenotypes(omim_id)
            if result and result.get("phenotypes"):
                return result
        # fallback to Monarch
        result = _monarch_disease_phenotypes(uid)
        if result:
            return result
        return {"disease_id": uid, "disease_name": "", "phenotypes": [], "genes": []}

    # OMIM / ORPHA — direct HPO API
    result = _hpo_api_disease_phenotypes(uid)
    return result or {"disease_id": uid, "disease_name": "", "phenotypes": [], "genes": []}


# ── HPO term → genes ───────────────────────────────────────────────────────

def _hpo_term_genes_hpo_api(hpo_id: str) -> list[str]:
    try:
        data = _get(f"{HPO_BASE}/term/{hpo_id}/genes", params={"max": 300})
        genes = data.get("genes") or data.get("geneAssoc") or []
        return [
            g.get("geneSymbol") or g.get("symbol") or g.get("gene", "")
            for g in genes
            if g.get("geneSymbol") or g.get("symbol") or g.get("gene")
        ]
    except Exception:
        return []


def _hpo_term_genes_monarch(hpo_id: str) -> list[str]:
    try:
        data = _get(f"{MONARCH_BASE}/association/phenotype/{hpo_id}/gene",
                    params={"limit": 300})
        items = data.get("items") or data.get("associations") or []
        syms = []
        for item in items:
            subj = item.get("subject") or item.get("gene") or {}
            sym = subj.get("symbol") or subj.get("label") or subj.get("name", "")
            if sym:
                syms.append(sym)
        return syms
    except Exception:
        return []


def get_term_genes(hpo_id: str) -> list[str]:
    """Get gene symbols associated with a specific HPO term (HP:xxxxxxx)."""
    syms = _hpo_term_genes_hpo_api(hpo_id)
    if not syms:
        syms = _hpo_term_genes_monarch(hpo_id)
    return syms


# ── Main evaluation function ───────────────────────────────────────────────

def evaluate_ppi_hpo_overlap(
    gene: str,
    ppi_partners: list[str],
    disease_name: str,
    omim_id: str = None,
    mondo_id: str = None,
    max_phenotypes: int = 30,
) -> dict:
    """Evaluate overlap between PPI partners and HPO symptom-associated genes.

    Args:
        gene: target gene symbol
        ppi_partners: list of PPI partner gene symbols (hub-filtered)
        disease_name: disease name (fallback if no ID given)
        omim_id: OMIM ID string like "OMIM:201910"
        mondo_id: MONDO ID like "MONDO_0008728" or "MONDO:0008728"
        max_phenotypes: maximum HPO terms to query

    Returns dict with per_term overlap and summary statistics.
    """
    ppi_set = {p.upper() for p in ppi_partners if p}

    # Resolve disease ID
    disease_id = omim_id or mondo_id
    disease_label = disease_name

    # Get phenotypes
    try:
        disease_data = get_disease_phenotypes(disease_id) if disease_id else \
                       {"disease_id": "", "disease_name": "", "phenotypes": [], "genes": []}
    except Exception as e:
        return {"error": f"HPO API error: {e}", "disease_id": disease_id or ""}

    phenotypes = disease_data.get("phenotypes", [])[:max_phenotypes]
    disease_label = disease_data.get("disease_name") or disease_label

    if not phenotypes:
        return {
            "error": "HPO症状データ取得失敗 (疾患IDが見つからないか症状なし)",
            "disease_id": disease_id or "",
            "disease_name": disease_label,
        }

    # Fetch gene sets per term in parallel
    def _fetch_term(pheno):
        hpo_id = pheno["hpo_id"]
        genes  = get_term_genes(hpo_id)
        return hpo_id, genes

    term_genes: dict[str, list[str]] = {}
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(_fetch_term, p): p for p in phenotypes}
        for fut in as_completed(futs):
            try:
                hpo_id, genes = fut.result()
                term_genes[hpo_id] = genes
            except Exception:
                pass

    # Per-term overlap
    per_term = []
    all_hpo_genes: set[str] = set()
    for pheno in phenotypes:
        hpo_id = pheno["hpo_id"]
        tgenes = term_genes.get(hpo_id, [])
        tgenes_upper = {g.upper() for g in tgenes if g}
        all_hpo_genes.update(tgenes_upper)
        overlap = sorted(ppi_set & tgenes_upper)
        per_term.append({
            "hpo_id":         hpo_id,
            "name":           pheno["name"],
            "frequency":      pheno.get("frequency", ""),
            "hpo_gene_count": len(tgenes_upper),
            "overlap_genes":  overlap,
            "overlap_count":  len(overlap),
        })

    per_term.sort(key=lambda x: x["overlap_count"], reverse=True)

    # Aggregate
    overlap_genes = ppi_set & all_hpo_genes
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

    return {
        "disease_id":     disease_id or "",
        "disease_name":   disease_label,
        "hpo_term_count": len(phenotypes),
        "per_term":       per_term,
        "summary": {
            "total_hpo_genes":   total_hpo_genes,
            "ppi_partner_count": len(ppi_set),
            "overlap_genes":     sorted(list(overlap_genes)),
            "overlap_count":     len(overlap_genes),
            "overlap_score":     overlap_score,
            "target_in_hpo":     gene.upper() in all_hpo_genes,
            "top_genes":         top_genes,
        },
    }
