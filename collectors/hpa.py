"""Human Protein Atlas — tissue/cell/pathology expression (CC BY-SA 4.0, 商用利用可).

タンパク質レベルの発現・病理データ。mRNA(GTEx)との比較に有用。
API: https://www.proteinatlas.org/{gene}.json
"""
import requests

HPA_BASE = "https://www.proteinatlas.org"


def get_expression_profile(gene_symbol: str) -> dict:
    """Return HPA tissue expression and disease pathology data.

    Returns:
        {
          "tissue_expression": [{tissue, level, reliability}],
          "subcellular":       [str],
          "cancer_expression": [{cancer, high_pct, low_pct}],
          "secretome_class":   str,
          "protein_class":     [str],
          "url":               str,
        }
    """
    url = f"{HPA_BASE}/{gene_symbol}.json"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 404:
            return {"error": f"{gene_symbol} not found in HPA"}
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        return {"error": str(e)}

    # 組織発現（タンパク質レベル）
    tissue_expr = []
    for entry in (data.get("Normal tissue") or []):
        level = entry.get("Level", "")
        if level in ("High", "Medium", "Low"):
            tissue_expr.append({
                "tissue":      entry.get("Tissue", ""),
                "cell_type":   entry.get("Cell type", ""),
                "level":       level,
                "reliability": entry.get("Reliability", ""),
            })

    # 細胞内局在
    subcellular = list({
        loc.get("Location", "")
        for loc in (data.get("Subcellular location") or [])
        if loc.get("Location")
    })

    # がん病理発現
    cancer_expr = []
    for entry in (data.get("Pathology") or [])[:10]:
        cancer_expr.append({
            "cancer":    entry.get("Cancer", ""),
            "high_pct":  entry.get("High", 0),
            "medium_pct":entry.get("Medium", 0),
            "low_pct":   entry.get("Low", 0),
            "not_det":   entry.get("Not detected", 0),
        })

    protein_class  = data.get("Protein class", [])
    secretome      = data.get("Secretome location", "")

    # 高発現組織のみ絞り込み（上位10件）
    high_tissues = [t for t in tissue_expr if t["level"] == "High"][:10]
    if not high_tissues:
        high_tissues = tissue_expr[:10]

    return {
        "tissue_expression": high_tissues,
        "subcellular":       subcellular,
        "cancer_expression": cancer_expr,
        "protein_class":     protein_class if isinstance(protein_class, list) else [protein_class],
        "secretome_class":   secretome,
        "url": f"{HPA_BASE}/{gene_symbol}",
    }
