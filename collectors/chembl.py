"""ChEMBL REST API (CC BY-SA 3.0) — drug and bioactivity data."""
import json
import time
import requests
from pathlib import Path

BASE = "https://www.ebi.ac.uk/chembl/api/data"

# ── キャッシュ設定 ─────────────────────────────────────────────────────────────
_CACHE_DIR = Path(__file__).parent.parent / "ppi_cache" / "chembl"
_CACHE_TTL = 7 * 24 * 3600  # 7日間

def _cache_path(gene_symbol: str) -> Path:
    return _CACHE_DIR / f"{gene_symbol.upper()}.json"

def _load_cache(gene_symbol: str):
    p = _cache_path(gene_symbol)
    if p.exists() and time.time() - p.stat().st_mtime < _CACHE_TTL:
        return json.loads(p.read_text(encoding="utf-8"))
    return None

def _save_cache(gene_symbol: str, data):
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    _cache_path(gene_symbol).write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")


def _find_target_chembl_id(gene_symbol: str) -> str | None:
    """遺伝子シンボルから ChEMBL target ID を取得する。

    1st: /target/search（全文検索）
    2nd: /target?target_synonym__icontains=（フォールバック）
    3rd: /target_component?accession=（UniProt 経由）
    """
    # --- 方法1: 全文検索（SINGLE PROTEIN・遺伝子名完全一致を優先） ---
    try:
        r = requests.get(f"{BASE}/target/search", params={
            "q": gene_symbol, "format": "json", "limit": 15
        }, timeout=10)
        if r.status_code == 200:
            human = [t for t in r.json().get("targets", [])
                     if t.get("organism") == "Homo sapiens"]
            single = [t for t in human if t.get("target_type") == "SINGLE PROTEIN"]
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

    # --- 方法2: synonym フィルター ---
    try:
        r2 = requests.get(f"{BASE}/target", params={
            "target_synonym__icontains": gene_symbol,
            "organism": "Homo sapiens",
            "target_type": "SINGLE PROTEIN",
            "format": "json", "limit": 5,
        }, timeout=10)
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
            }, timeout=10)
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
    """Return approved/clinical-stage drugs that target the given gene."""
    cached = _load_cache(gene_symbol)
    if cached is not None:
        return cached

    chembl_id = _find_target_chembl_id(gene_symbol)
    if not chembl_id:
        _save_cache(gene_symbol, [])
        return []

    # mechanism テーブルから承認薬を取得
    try:
        r2 = requests.get(f"{BASE}/mechanism", params={
            "target_chembl_id": chembl_id,
            "format": "json",
            "limit": max_results,
        }, timeout=15)
        r2.raise_for_status()
        mechanisms = r2.json().get("mechanisms", [])
    except Exception:
        mechanisms = []

    mech_by_mol = {}
    for mech in mechanisms:
        mol_id = mech.get("molecule_chembl_id")
        if mol_id and mol_id not in mech_by_mol:
            mech_by_mol[mol_id] = mech

    drugs = []
    if mech_by_mol:
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
                }, timeout=15)
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

    # 承認薬がない場合は臨床フェーズ化合物にフォールバック（重いので件数を絞る）
    if not drugs:
        drugs = _clinical_candidates(chembl_id, max_n=min(max_results, 30))

    result = drugs[:max_results]
    _save_cache(gene_symbol, result)
    return result


def _clinical_candidates(target_chembl_id: str, max_n: int = 30) -> list[dict]:
    """mechanism テーブルに承認薬が無いターゲット向け: 臨床フェーズ化合物を返す。"""
    try:
        r = requests.get(f"{BASE}/activity", params={
            "target_chembl_id": target_chembl_id,
            "pchembl_value__isnull": "false",
            "format": "json", "limit": max_n * 2,  # 絞り込み後に max_n 件確保
        }, timeout=20)
        r.raise_for_status()
        acts = r.json().get("activities", [])
    except Exception:
        return []

    mol_ids = list({a.get("molecule_chembl_id") for a in acts if a.get("molecule_chembl_id")})[:max_n * 2]

    mol_by_id = {}
    CHUNK = 50
    for i in range(0, len(mol_ids), CHUNK):
        chunk = mol_ids[i:i + CHUNK]
        try:
            r3 = requests.get(f"{BASE}/molecule", params={
                "molecule_chembl_id__in": ",".join(chunk),
                "format": "json",
                "limit": len(chunk),
            }, timeout=15)
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
        }, timeout=10)
        if r.status_code != 200:
            return {"has_safety_flag": False}
        targets = r.json().get("targets", [])
    except Exception:
        return {"has_safety_flag": False}

    for t in targets:
        for comp in (t.get("target_components") or []):
            for prop in (comp.get("target_component_synonyms") or []):
                if "safety" in (prop.get("component_synonym") or "").lower():
                    return {"has_safety_flag": True, "details": prop}
    return {"has_safety_flag": False}
