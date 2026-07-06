"""Human Protein Atlas — tissue/cell expression (CC BY-SA 4.0, 商用利用可).

API: https://www.proteinatlas.org/{ENSG_ID}.json
"""
import requests
from collectors._ensembl import resolve_ensg

HPA_BASE = "https://www.proteinatlas.org"


def get_expression_profile(gene_symbol: str) -> dict:
    """Return HPA tissue expression and subcellular localisation.

    Returns:
        {
          "tissue_expression": [{tissue, level}],   # RNA nTPM, sorted descending
          "protein_tissue":    [{tissue, level}],   # protein intensity (qualitative)
          "subcellular":       [str],
          "protein_class":     [str],
          "url":               str,
        }
    """
    ensg = resolve_ensg(gene_symbol)
    if not ensg:
        return {"error": f"Could not resolve ENSG for {gene_symbol}"}

    url = f"{HPA_BASE}/{ensg}.json"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 404:
            return {"error": f"{gene_symbol} ({ensg}) not found in HPA"}
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        return {"error": str(e)}

    # RNA tissue expression — nTPM dict {tissue: str_or_float}
    rna_tissue = data.get("RNA tissue specific nTPM") or {}
    tissue_expr = []
    for tissue, val in rna_tissue.items():
        try:
            tpm = float(val)
        except (TypeError, ValueError):
            tpm = 0.0
        tissue_expr.append({"tissue": tissue, "level": f"{tpm:.1f} nTPM", "tpm": tpm})
    tissue_expr.sort(key=lambda x: x["tpm"], reverse=True)
    tissue_expr = tissue_expr[:15]

    # Protein tissue intensity (qualitative)
    prot_tissue = data.get("Protein tissue specific Intensity") or {}
    prot_expr = [{"tissue": t, "level": v} for t, v in prot_tissue.items()]

    # Subcellular localisation
    subcellular = []
    for k in ("Subcellular main location", "Subcellular additional location", "Subcellular location"):
        val = data.get(k)
        if isinstance(val, list):
            subcellular.extend(val)
        elif isinstance(val, str) and val:
            subcellular.append(val)
    subcellular = list(dict.fromkeys(subcellular))

    protein_class = data.get("Protein class") or []
    if isinstance(protein_class, str):
        protein_class = [protein_class]

    return {
        "tissue_expression": tissue_expr,
        "protein_tissue":    prot_expr,
        "subcellular":       subcellular,
        "protein_class":     protein_class,
        "url": f"{HPA_BASE}/{ensg}",
    }
