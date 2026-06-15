"""ChEMBL REST API (CC BY-SA 3.0) — drug and bioactivity data."""
import requests

BASE = "https://www.ebi.ac.uk/chembl/api/data"

def _find_target_chembl_id(gene_symbol: str) -> str | None:
    """遺伝子シンボルから ChEMBL target ID を取得する。

    1st: /target/search（全文検索）
    2nd: /target?pref_name__icontains=（500 時のフォールバック）
    3rd: /target_component?accession=（UniProt 経由）
    """
    # --- 方法1: 全文検索 ---
    try:
        r = requests.get(f"{BASE}/target/search", params={
            "q": gene_symbol, "format": "json", "limit": 5
        }, timeout=20)
        if r.status_code == 200:
            human = [t for t in r.json().get("targets", [])
                     if t.get("organism") == "Homo sapiens"]
            if human:
                return human[0].get("target_chembl_id")
    except Exception:
        pass

    # --- 方法2: gene_name フィルター ---
    try:
        r2 = requests.get(f"{BASE}/target", params={
            "target_synonym__icontains": gene_symbol,
            "organism": "Homo sapiens",
            "target_type": "SINGLE PROTEIN",
            "format": "json", "limit": 5,
        }, timeout=20)
        if r2.status_code == 200:
            targets = r2.json().get("targets", [])
            if targets:
                return targets[0].get("target_chembl_id")
    except Exception:
        pass

    # --- 方法3: UniProt accession 経由 ---
    try:
        from collectors.uniprot import get_protein_info
        info = get_protein_info(gene_symbol)
        uniprot_id = info.get("uniprot_id", "")
        if uniprot_id:
            r3 = requests.get(f"{BASE}/target_component", params={
                "accession": uniprot_id, "format": "json", "limit": 3,
            }, timeout=20)
            if r3.status_code == 200:
                comps = r3.json().get("target_components", [])
                for comp in comps:
                    for link in (comp.get("targets") or []):
                        tid = link.get("target_chembl_id")
                        if tid:
                            return tid
    except Exception:
        pass

    return None


def get_drugs_for_target(gene_symbol: str) -> list[dict]:
    """Return approved/clinical-stage drugs that target the given gene."""
    chembl_id = _find_target_chembl_id(gene_symbol)
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
    try:
        r = requests.get(f"{BASE}/target/search", params={
            "q": gene_symbol, "format": "json", "limit": 3
        }, timeout=20)
        if r.status_code != 200:
            return {"has_safety_flag": False}
        targets = r.json().get("targets", [])
    except Exception:
        return {"has_safety_flag": False}

    targets = targets  # 以下の既存ロジックに続く
    for t in targets:
        props = t.get("target_components", [])
        for comp in props:
            for prop in comp.get("target_component_synonyms", []):
                if "safety" in prop.get("component_synonym", "").lower():
                    return {"has_safety_flag": True, "details": prop}
    return {"has_safety_flag": False}
