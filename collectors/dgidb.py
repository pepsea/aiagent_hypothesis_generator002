"""DGIdb — Drug-Gene Interaction Database (Apache 2.0, 商用利用可).

ChEMBL を補完する薬剤-遺伝子相互作用データ。リポジショニング候補発掘に有用。
GraphQL API: https://dgidb.org/api/graphql
"""
import requests

DGIDB_API = "https://dgidb.org/api/graphql"

INTERACTIONS_QUERY = """
query($gene: String!) {
  genes(names: [$gene]) {
    nodes {
      name
      interactions {
        drug {
          name
          approved
          drugAttributes { name value }
        }
        interactionScore
        interactionTypes { type directionality }
        sources { sourceDbName fullName }
        pmids
      }
    }
  }
}
"""


def get_interactions(gene_symbol: str, max_results: int = 20) -> list[dict]:
    """Return drug-gene interactions from DGIdb.

    Returns:
        [{drug_name, approved, interaction_type, directionality,
          score, sources, pmids}]
    """
    try:
        r = requests.post(
            DGIDB_API,
            json={"query": INTERACTIONS_QUERY, "variables": {"gene": gene_symbol}},
            timeout=20,
        )
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        return []

    nodes = (data.get("data", {}).get("genes", {}).get("nodes") or [])
    if not nodes:
        return []

    interactions = nodes[0].get("interactions") or []
    results = []
    seen = set()

    for ix in interactions:
        drug = ix.get("drug") or {}
        name = (drug.get("name") or "").upper()
        if not name or name in seen:
            continue
        seen.add(name)

        i_types = ix.get("interactionTypes") or []
        i_type  = i_types[0].get("type", "") if i_types else ""
        directionality = i_types[0].get("directionality", "") if i_types else ""

        sources = [s.get("sourceDbName", "") for s in (ix.get("sources") or [])]
        pmids   = (ix.get("pmids") or [])[:3]

        results.append({
            "drug_name":      drug.get("name", ""),
            "approved":       drug.get("approved", False),
            "interaction_type": i_type,
            "directionality": directionality,
            "score":          ix.get("interactionScore"),
            "sources":        sources,
            "pmids":          pmids,
        })

    # approved 優先、score 降順
    results.sort(key=lambda x: (x["approved"] is True, x["score"] or 0), reverse=True)
    return results[:max_results]
