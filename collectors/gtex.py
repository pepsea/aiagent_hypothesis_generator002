"""GTEx — tissue gene expression (CC BY 4.0, 商用利用可).

組織別発現プロファイルからモダリティ選択・副作用リスクを推定。
API: https://gtexportal.org/api/v2
"""
import requests

GTEX_API = "https://gtexportal.org/api/v2"

# 疾患関連・安全性評価に重要な組織リスト
KEY_TISSUES = {
    "Heart_Left_Ventricle":    "Heart (left ventricle)",
    "Liver":                   "Liver",
    "Kidney_Cortex":           "Kidney (cortex)",
    "Brain_Frontal_Cortex_Ba9":"Brain (frontal cortex)",
    "Lung":                    "Lung",
    "Muscle_Skeletal":         "Skeletal muscle",
    "Colon_Sigmoid":           "Colon",
    "Whole_Blood":             "Whole blood",
    "Adipose_Subcutaneous":    "Adipose tissue",
    "Skin_Sun_Exposed_Lower_leg": "Skin",
}


def get_tissue_expression(gene_symbol: str, top_n: int = 10) -> dict:
    """Return median TPM expression per tissue for the gene.

    Returns:
        {
          "top_tissues":   [{tissue, tpm}],   # 全組織中上位 top_n
          "key_tissues":   [{tissue, tpm}],   # 安全性関連組織
          "max_tissue":    str,
          "max_tpm":       float,
          "url": str,
        }
    """
    # Ensembl ID → GTEx gene ID 変換 (遺伝子シンボルで検索)
    try:
        r = requests.get(f"{GTEX_API}/reference/gene", params={
            "geneSymbol": gene_symbol,
            "gencodeVersion": "v26",
            "genomeBuild": "GRCh38/hg38",
        }, timeout=15)
        r.raise_for_status()
        genes = r.json().get("data", [])
        if not genes:
            return {"error": f"Gene {gene_symbol} not found in GTEx"}
        versioned_id = genes[0].get("gencodeId", "")
    except Exception as e:
        return {"error": str(e)}

    # 組織別発現取得
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

    # 全組織ソート
    all_tissues = sorted(
        [{"tissue": r.get("tissueSiteDetailId", ""), "tpm": r.get("median", 0)} for r in records],
        key=lambda x: x["tpm"], reverse=True,
    )

    top_tissues = all_tissues[:top_n]

    # 安全性関連組織の発現
    tpm_map = {r.get("tissueSiteDetailId", ""): r.get("median", 0) for r in records}
    key_tissues = [
        {"tissue": label, "tpm": tpm_map.get(tid, 0)}
        for tid, label in KEY_TISSUES.items()
    ]
    key_tissues.sort(key=lambda x: x["tpm"], reverse=True)

    max_t = top_tissues[0] if top_tissues else {}

    return {
        "top_tissues": top_tissues,
        "key_tissues": key_tissues,
        "max_tissue":  max_t.get("tissue", ""),
        "max_tpm":     max_t.get("tpm", 0),
        "url": f"https://gtexportal.org/home/gene/{gene_symbol}",
    }
