"""IntAct protein interaction database (EMBL-EBI, CC BY 4.0).

String代替として使用。商用利用可。
"""
import requests

BASE = "https://www.ebi.ac.uk/intact/ws/interaction"

def get_interactions(gene_symbol: str, species: int = 9606, max_results: int = 20) -> list[dict]:
    """Return top interactors for a gene from IntAct."""
    r = requests.get(f"{BASE}/findInteractions/{gene_symbol}", params={
        "page": 0, "pageSize": max_results,
        "query": f"species:{species}",
    }, timeout=20)

    if r.status_code == 404:
        return []
    r.raise_for_status()

    data = r.json()
    interactions = []
    for item in data.get("content", []):
        participants = item.get("participants", [])

        # participants はdictのリストまたは文字列のリストの場合がある
        names = []
        for p in participants:
            if isinstance(p, dict):
                alias = p.get("preferredName", "") or p.get("interactorAc", "")
            else:
                alias = str(p)
            if alias:
                names.append(alias)

        partners = [n for n in names if gene_symbol.upper() not in n.upper()]

        pubs = item.get("publications", [])
        pubmed_ids = []
        for p in pubs:
            if isinstance(p, dict):
                pubmed_ids.append(p.get("pubmedId", ""))
            else:
                pubmed_ids.append(str(p))

        interactions.append({
            "interaction_id": item.get("interactionAc", ""),
            "partners": partners,
            "detection_method": (item.get("detectionMethod") or {}).get("shortName", "") if isinstance(item.get("detectionMethod"), dict) else "",
            "interaction_type": (item.get("interactionType") or {}).get("shortName", "") if isinstance(item.get("interactionType"), dict) else "",
            "confidence": item.get("intactScore", None),
            "pubmed_ids": pubmed_ids,
        })

    return interactions


def get_top_interactors(gene_symbol: str, top_n: int = 10) -> list[str]:
    """Return list of top interactor gene names."""
    interactions = get_interactions(gene_symbol, max_results=50)
    partners = []
    for ix in interactions:
        partners.extend(ix["partners"])
    # Count frequency
    from collections import Counter
    counts = Counter(partners)
    return [name for name, _ in counts.most_common(top_n)]
