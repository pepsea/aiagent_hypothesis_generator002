"""AlphaFold DB — predicted protein structures (CC BY 4.0, 商用利用可).

ドラッガビリティ評価（pocket有無）・モダリティ選択の根拠。
pLDDT スコアで構造信頼性を定量化。
API: https://alphafold.ebi.ac.uk/api
"""
import requests

AF_API = "https://alphafold.ebi.ac.uk/api"


def get_structure_info(gene_symbol: str, uniprot_id: str = "") -> dict:
    """Return AlphaFold structure confidence and links for the protein.

    Tries UniProt ID directly; falls back to UniProt search by gene symbol.

    Returns:
        {
          "uniprot_id":     str,
          "entry_id":       str,   # e.g. AF-P05067-F1
          "mean_plddt":     float, # 0-100; >90=very high, 70-90=high, 50-70=low, <50=very low
          "confidence":     str,   # interpretation
          "pdb_url":        str,
          "view_url":       str,
        }
    """
    uid = uniprot_id

    # UniProt ID が未指定の場合は UniProt REST で解決
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
        return {"error": f"UniProt ID not found for {gene_symbol}"}

    # AlphaFold エントリ取得
    try:
        r = requests.get(f"{AF_API}/prediction/{uid}", timeout=15)
        if r.status_code == 404:
            return {"error": f"AlphaFold entry not found for {uid}"}
        r.raise_for_status()
        entries = r.json()
        if not entries:
            return {"error": "No AlphaFold entry"}
        entry = entries[0]
    except Exception as e:
        return {"error": str(e)}

    # API field renamed: meanPlddt → globalMetricValue
    mean_plddt = entry.get("globalMetricValue")
    if mean_plddt is None:
        mean_plddt = entry.get("meanPlddt")
    entry_id   = entry.get("entryId", "")
    pdb_url    = entry.get("pdbUrl", "")
    view_url   = f"https://alphafold.ebi.ac.uk/entry/{uid}"

    # pLDDT 解釈
    if mean_plddt is not None:
        if mean_plddt >= 90:
            confidence = "Very high (≥90): reliable backbone, suitable for structure-based drug design"
        elif mean_plddt >= 70:
            confidence = "High (70–90): mostly reliable, moderate confidence"
        elif mean_plddt >= 50:
            confidence = "Low (50–70): disordered regions likely, limited SBDD utility"
        else:
            confidence = "Very low (<50): intrinsically disordered, consider non-SBDD modalities"
    else:
        confidence = "Unknown"

    return {
        "uniprot_id": uid,
        "entry_id":   entry_id,
        "mean_plddt": mean_plddt,
        "confidence": confidence,
        "pdb_url":    pdb_url,
        "view_url":   view_url,
    }
