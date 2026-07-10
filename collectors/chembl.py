"""ChEMBL REST API (CC BY-SA 3.0) — drug and bioactivity data."""
import requests

BASE = "https://www.ebi.ac.uk/chembl/api/data"

def _find_target_chembl_id(gene_symbol: str) -> str | None:
    """遺伝子シンボルから ChEMBL target ID を取得する。

    1st: /target/search（全文検索）
    2nd: /target?pref_name__icontains=（500 時のフォールバック）
    3rd: /target_component?accession=（UniProt 経由）
    """
    # --- 方法1: 全文検索（SINGLE PROTEIN・遺伝子名完全一致を優先） ---
    try:
        r = requests.get(f"{BASE}/target/search", params={
            "q": gene_symbol, "format": "json", "limit": 15
        }, timeout=20)
        if r.status_code == 200:
            human = [t for t in r.json().get("targets", [])
                     if t.get("organism") == "Homo sapiens"]
            # 単一タンパク質を最優先
            single = [t for t in human if t.get("target_type") == "SINGLE PROTEIN"]
            # 遺伝子シンボルがコンポーネント synonym に完全一致するもの
            def _matches_symbol(t):
                for comp in (t.get("target_components") or []):
                    for syn in (comp.get("target_component_synonyms") or []):
                        if (syn.get("component_synonym") or "").upper() == gene_symbol.upper():
                            return True
                return False
            exact = [t for t in single if _matches_symbol(t)]
            for pool in (exact, single, human):
                if pool:
                    return pool[0].get("target_chembl_id")
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


def get_drugs_for_target(gene_symbol: str, max_results: int = 100) -> list[dict]:
    """Return approved/clinical-stage drugs that target the given gene (最大 max_results 件)。"""
    chembl_id = _find_target_chembl_id(gene_symbol)
    if not chembl_id:
        return []

    # Get approved drugs for this target
    r2 = requests.get(f"{BASE}/mechanism", params={
        "target_chembl_id": chembl_id,
        "format": "json",
        "limit": max_results,
    }, timeout=20)
    r2.raise_for_status()

    mechanisms = r2.json().get("mechanisms", [])
    mech_by_mol = {}
    for mech in mechanisms:
        mol_id = mech.get("molecule_chembl_id")
        if mol_id and mol_id not in mech_by_mol:
            mech_by_mol[mol_id] = mech

    drugs = []
    if mech_by_mol:
        # molecule 詳細を1件ずつ取得すると件数が増えるほどN+1で遅くなるため、
        # molecule_chembl_id__in でバッチ取得する（ChEMBL API がサポート）
        mol_ids = list(mech_by_mol.keys())
        mol_by_id = {}
        CHUNK = 50
        for i in range(0, len(mol_ids), CHUNK):
            chunk = mol_ids[i:i + CHUNK]
            try:
                r3 = requests.get(f"{BASE}/molecule", params={
                    "molecule_chembl_id__in": ",".join(chunk),
                    "format": "json",
                    "limit": len(chunk),
                }, timeout=20)
                r3.raise_for_status()
                for mol in r3.json().get("molecules", []):
                    mol_by_id[mol.get("molecule_chembl_id")] = mol
            except Exception:
                pass

        for mol_id, mech in mech_by_mol.items():
            mol = mol_by_id.get(mol_id)
            if mol:
                drugs.append({
                    "chembl_id": mol_id,
                    "name": mol.get("pref_name", mol_id),
                    "max_phase": mol.get("max_phase"),
                    "molecule_type": mol.get("molecule_type", ""),
                    "mechanism": mech.get("mechanism_of_action", ""),
                    "action_type": mech.get("action_type", ""),
                    "indication": mech.get("disease_efficacy", False),
                })
            else:
                drugs.append({
                    "chembl_id": mol_id,
                    "name": mol_id,
                    "mechanism": mech.get("mechanism_of_action", ""),
                    "action_type": mech.get("action_type", ""),
                })

    # 承認薬 mechanism が無い場合、活性を持つ臨床フェーズ化合物にフォールバック
    if not drugs:
        drugs = _clinical_candidates(chembl_id, max_n=max_results)

    return drugs[:max_results]


def _clinical_candidates(target_chembl_id: str, max_n: int = 100) -> list[dict]:
    """mechanism テーブルに承認薬が無いターゲット向け:
    活性データを持つ臨床フェーズ (max_phase>=1) の名前付き分子を返す。"""
    try:
        r = requests.get(f"{BASE}/activity", params={
            "target_chembl_id": target_chembl_id,
            "pchembl_value__isnull": "false",
            "format": "json", "limit": max(max_n, 100),
        }, timeout=30)
        r.raise_for_status()
        acts = r.json().get("activities", [])
    except Exception:
        return []

    mol_ids = list({a.get("molecule_chembl_id") for a in acts if a.get("molecule_chembl_id")})

    # molecule 詳細をバッチ取得（1件ずつだと max_n=100 で N+1 が重くなるため）
    mol_by_id = {}
    CHUNK = 50
    for i in range(0, len(mol_ids), CHUNK):
        chunk = mol_ids[i:i + CHUNK]
        try:
            r3 = requests.get(f"{BASE}/molecule", params={
                "molecule_chembl_id__in": ",".join(chunk),
                "format": "json",
                "limit": len(chunk),
            }, timeout=20)
            r3.raise_for_status()
            for mol in r3.json().get("molecules", []):
                mol_by_id[mol.get("molecule_chembl_id")] = mol
        except Exception:
            pass

    drugs, seen = [], set()
    for mol_id in mol_ids:
        mol = mol_by_id.get(mol_id)
        if not mol:
            continue
        phase = mol.get("max_phase")
        name = mol.get("pref_name")
        # 臨床フェーズかつ名前付きのもののみ
        if not name or phase in (None, 0, "0") or name in seen:
            continue
        seen.add(name)
        drugs.append({
            "chembl_id": mol_id,
            "name": name,
            "max_phase": phase,
            "molecule_type": mol.get("molecule_type", ""),
            "mechanism": "(bioactivity; not an approved indication)",
            "action_type": "",
            "indication": False,
        })

    drugs.sort(key=lambda d: d.get("max_phase") or 0, reverse=True)
    return drugs[:max_n]


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
