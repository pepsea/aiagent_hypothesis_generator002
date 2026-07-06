"""GTEx — tissue gene expression (CC BY 4.0, 商用利用可).

API: https://gtexportal.org/api/v2  (v2 requires Ensembl geneId)
"""
import requests
from collectors._ensembl import resolve_ensg

GTEX_API = "https://gtexportal.org/api/v2"

KEY_TISSUES = {
    "Heart_Left_Ventricle":     "Heart (left ventricle)",
    "Liver":                    "Liver",
    "Kidney_Cortex":            "Kidney (cortex)",
    "Brain_Frontal_Cortex_Ba9": "Brain (frontal cortex)",
    "Brain_Substantia_nigra":   "Brain (substantia nigra)",
    "Lung":                     "Lung",
    "Muscle_Skeletal":          "Skeletal muscle",
    "Colon_Sigmoid":            "Colon",
    "Whole_Blood":              "Whole blood",
    "Adipose_Subcutaneous":     "Adipose tissue",
}


def get_tissue_expression(gene_symbol: str, top_n: int = 10) -> dict:
    """Return median TPM expression per tissue for the gene."""
    ensembl_id = resolve_ensg(gene_symbol)
    if not ensembl_id:
        return {"error": f"Could not resolve Ensembl ID for {gene_symbol}"}

    # Step 1: versioned gencodeId
    try:
        r = requests.get(f"{GTEX_API}/reference/gene",
                         params={"geneId": ensembl_id}, timeout=15)
        r.raise_for_status()
        genes = r.json().get("data", [])
        if not genes:
            return {"error": f"{gene_symbol} not found in GTEx"}
        versioned_id = genes[0].get("gencodeId", "")
    except Exception as e:
        return {"error": str(e)}

    # Step 2: median expression per tissue
    try:
        r2 = requests.get(f"{GTEX_API}/expression/medianGeneExpression", params={
            "gencodeId": versioned_id,
            "datasetId": "gtex_v8",
        }, timeout=20)
        r2.raise_for_status()
        records = r2.json().get("data", [])
    except Exception as e:
        return {"error": str(e)}

    if not records:
        return {"error": "No expression data"}

    all_tissues = sorted(
        [{"tissue": rec.get("tissueSiteDetailId", ""), "tpm": rec.get("median", 0)}
         for rec in records],
        key=lambda x: x["tpm"], reverse=True,
    )
    top_tissues = all_tissues[:top_n]
    tpm_map = {rec.get("tissueSiteDetailId", ""): rec.get("median", 0) for rec in records}
    key_tissues = sorted(
        [{"tissue": label, "tpm": tpm_map.get(tid, 0)} for tid, label in KEY_TISSUES.items()],
        key=lambda x: x["tpm"], reverse=True,
    )
    max_t = top_tissues[0] if top_tissues else {}

    return {
        "top_tissues": top_tissues,
        "key_tissues": key_tissues,
        "max_tissue":  max_t.get("tissue", ""),
        "max_tpm":     max_t.get("tpm", 0),
        "url": f"https://gtexportal.org/home/gene/{gene_symbol}",
    }
