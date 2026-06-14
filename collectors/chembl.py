"""ChEMBL REST API (CC BY-SA 3.0) — drug and bioactivity data."""
import requests

BASE = "https://www.ebi.ac.uk/chembl/api/data"

def get_drugs_for_target(gene_symbol: str) -> list[dict]:
    """Return approved/clinical-stage drugs that target the given gene."""
    # Search target by gene name
    r = requests.get(f"{BASE}/target/search", params={
        "q": gene_symbol, "format": "json", "limit": 5
    }, timeout=20)
    r.raise_for_status()

    targets = r.json().get("targets", [])
    human_targets = [t for t in targets if t.get("organism") == "Homo sapiens"]
    if not human_targets:
        return []

    chembl_id = human_targets[0].get("target_chembl_id")
    if not chembl_id:
        return []

    # Get approved drugs for this target
    r2 = requests.get(f"{BASE}/mechanism", params={
        "target_chembl_id": chembl_id,
        "format": "json",
        "limit": 20,
    }, timeout=20)
    r2.raise_for_status()

    mechanisms = r2.json().get("mechanisms", [])
    drugs = []
    seen = set()
    for mech in mechanisms:
        mol_id = mech.get("molecule_chembl_id")
        if mol_id in seen:
            continue
        seen.add(mol_id)

        # Get drug details
        try:
            r3 = requests.get(f"{BASE}/molecule/{mol_id}", params={"format": "json"}, timeout=10)
            r3.raise_for_status()
            mol = r3.json()
            drugs.append({
                "chembl_id": mol_id,
                "name": mol.get("pref_name", mol_id),
                "max_phase": mol.get("max_phase"),
                "molecule_type": mol.get("molecule_type", ""),
                "mechanism": mech.get("mechanism_of_action", ""),
                "action_type": mech.get("action_type", ""),
                "indication": mech.get("disease_efficacy", False),
            })
        except Exception:
            drugs.append({
                "chembl_id": mol_id,
                "name": mol_id,
                "mechanism": mech.get("mechanism_of_action", ""),
                "action_type": mech.get("action_type", ""),
            })

    return drugs


def get_toxicity_flags(gene_symbol: str) -> dict:
    """Check if target gene has toxicity-related flags in ChEMBL."""
    r = requests.get(f"{BASE}/target/search", params={
        "q": gene_symbol, "format": "json", "limit": 3
    }, timeout=20)
    r.raise_for_status()

    targets = r.json().get("targets", [])
    for t in targets:
        props = t.get("target_components", [])
        for comp in props:
            for prop in comp.get("target_component_synonyms", []):
                if "safety" in prop.get("component_synonym", "").lower():
                    return {"has_safety_flag": True, "details": prop}
    return {"has_safety_flag": False}
