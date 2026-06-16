"""gnomAD — population genetics & constraint scores (CC BY 4.0, 商用利用可).

制約スコア (pLI, LOEUF) でターゲットの必須性・毒性リスクを定量化。
GraphQL API: https://gnomad.broadinstitute.org/api
"""
import requests

GNOMAD_API = "https://gnomad.broadinstitute.org/api"

CONSTRAINT_QUERY = """
query GeneConstraint($geneSymbol: String!) {
  gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
    gene_id
    gnomad_constraint {
      exp_lof
      obs_lof
      lof_z
      pLI
      oe_lof
      oe_lof_upper
    }
    gnomad_missense_constraint: gnomad_constraint {
      exp_mis
      obs_mis
      mis_z
      oe_mis
      oe_mis_upper
    }
  }
}
"""


def get_constraint(gene_symbol: str) -> dict:
    """Return gnomAD constraint scores for the gene.

    Key metrics:
      pLI   : probability of loss-of-function intolerance (>0.9 = highly constrained)
      LOEUF : loss-of-function observed/expected upper bound fraction (<0.35 = constrained)
      oe_mis: missense observed/expected ratio
    """
    try:
        r = requests.post(
            GNOMAD_API,
            json={"query": CONSTRAINT_QUERY, "variables": {"geneSymbol": gene_symbol}},
            timeout=20,
        )
        r.raise_for_status()
        data = r.json().get("data", {}).get("gene") or {}
    except Exception as e:
        return {"error": str(e)}

    c = data.get("gnomad_constraint") or {}
    if not c:
        return {"error": "constraint data not found"}

    pli     = c.get("pLI")
    loeuf   = c.get("oe_lof_upper")
    lof_z   = c.get("lof_z")
    obs_lof = c.get("obs_lof")
    exp_lof = c.get("exp_lof")
    oe_mis  = c.get("oe_mis")

    # 必須性の解釈
    if pli is not None:
        if pli >= 0.9:
            essentiality = "High (pLI≥0.9: likely essential — elevated on-target toxicity risk)"
        elif pli >= 0.5:
            essentiality = "Moderate (pLI 0.5–0.9)"
        else:
            essentiality = "Low (pLI<0.5: tolerates LoF variants)"
    else:
        essentiality = "Unknown"

    return {
        "gene_id":      data.get("gene_id", ""),
        "pLI":          pli,
        "LOEUF":        loeuf,
        "lof_z":        lof_z,
        "obs_lof":      obs_lof,
        "exp_lof":      exp_lof,
        "oe_missense":  oe_mis,
        "essentiality": essentiality,
        "url": f"https://gnomad.broadinstitute.org/gene/{gene_symbol}?dataset=gnomad_r4",
    }
