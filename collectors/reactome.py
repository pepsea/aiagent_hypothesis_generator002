"""Reactome — pathway database (CC BY 4.0, 商用利用可).

因果関係付きパスウェイ情報。g:Profilerよりも詳細な経路コンテキスト。
ContentService API: https://reactome.org/ContentService
"""
import requests

REACTOME_API = "https://reactome.org/ContentService"


def get_pathways(gene_symbol: str, uniprot_id: str = "", top_n: int = 15) -> list[dict]:
    """Return top Reactome pathways for the gene (via UniProt accession).

    Returns:
        [{pathway_id, name, species, is_disease, url}]
    """
    uid = uniprot_id

    # UniProt ID 解決
    if not uid:
        try:
            r0 = requests.get("https://rest.uniprot.org/uniprotkb/search", params={
                "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
                "fields": "accession",
                "format": "json",
                "size": 1,
            }, timeout=10)
            r0.raise_for_status()
            results = r0.json().get("results", [])
            uid = results[0]["primaryAccession"] if results else ""
        except Exception:
            pass

    if not uid:
        return []

    # Reactome パスウェイ取得
    try:
        r = requests.get(
            f"{REACTOME_API}/data/pathways/low/entity/{uid}/allForms",
            params={"species": "9606"},
            timeout=20,
        )
        if r.status_code == 404:
            return []
        r.raise_for_status()
        pathways = r.json()
    except Exception:
        return []

    results = []
    for p in pathways[:top_n]:
        pid  = p.get("stId", "")
        name = p.get("displayName", "")
        results.append({
            "pathway_id": pid,
            "name":       name,
            "is_disease": "disease" in name.lower() or p.get("isInDisease", False),
            "url":        f"https://reactome.org/PathwayBrowser/#/{pid}",
        })

    return results


def get_disease_pathways(gene_symbol: str, uniprot_id: str = "") -> list[dict]:
    """Return only disease-annotated pathways."""
    all_pathways = get_pathways(gene_symbol, uniprot_id, top_n=50)
    return [p for p in all_pathways if p["is_disease"]]
