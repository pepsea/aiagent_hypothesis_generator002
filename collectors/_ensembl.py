"""Shared Ensembl gene symbol → ID resolver with simple in-process cache.

解決順:
  1. Ensembl REST API (rest.ensembl.org)
  2. MyGene.info API (フォールバック — Ensembl が到達不能な場合)
  3. OpenTargets GraphQL (最終フォールバック)
"""
import time
import requests

_CACHE: dict[str, str] = {}
ENSEMBL_REST = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
MYGENE_API   = "https://mygene.info/v3/query"
OT_API       = "https://api.platform.opentargets.org/api/v4/graphql"

_OT_SEARCH_Q = """
query($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits { id name }
  }
}
"""


def _resolve_via_mygene(gene_symbol: str) -> str:
    """MyGene.info で ENSG ID を解決する。"""
    try:
        r = requests.get(MYGENE_API, params={
            "q": f"symbol:{gene_symbol}",
            "species": "human",
            "fields": "ensembl.gene",
            "size": 1,
        }, timeout=10)
        r.raise_for_status()
        hits = r.json().get("hits", [])
        if not hits:
            return ""
        ensembl = hits[0].get("ensembl", {})
        # ensembl フィールドはリストの場合もある
        if isinstance(ensembl, list):
            ensembl = ensembl[0]
        return ensembl.get("gene", "")
    except Exception:
        return ""


def _resolve_via_opentargets(gene_symbol: str) -> str:
    """OpenTargets GraphQL で ENSG ID を解決する（最終手段）。"""
    try:
        r = requests.post(OT_API, json={
            "query": _OT_SEARCH_Q,
            "variables": {"q": gene_symbol},
        }, timeout=15)
        r.raise_for_status()
        hits = r.json().get("data", {}).get("search", {}).get("hits", [])
        for h in hits:
            if h.get("name", "").upper() == gene_symbol.upper():
                return h.get("id", "")
        return hits[0].get("id", "") if hits else ""
    except Exception:
        return ""


def resolve_ensg(gene_symbol: str, max_retries: int = 2) -> str:
    """Return Ensembl ENSG ID for a human gene symbol, or '' on failure."""
    key = gene_symbol.upper()
    if key in _CACHE:
        return _CACHE[key]

    # 1. Ensembl REST
    for attempt in range(max_retries):
        try:
            r = requests.get(
                f"{ENSEMBL_REST}/{gene_symbol}",
                headers={"Content-Type": "application/json"},
                timeout=10,
            )
            if r.status_code == 404:
                break  # 遺伝子が存在しない → フォールバックへ
            r.raise_for_status()
            ensg = r.json().get("id", "")
            if ensg:
                _CACHE[key] = ensg
                return ensg
        except Exception:
            if attempt < max_retries - 1:
                time.sleep(1)

    # 2. MyGene.info フォールバック
    ensg = _resolve_via_mygene(gene_symbol)
    if ensg:
        _CACHE[key] = ensg
        return ensg

    # 3. OpenTargets フォールバック（最も確実）
    ensg = _resolve_via_opentargets(gene_symbol)
    if ensg:
        _CACHE[key] = ensg
        return ensg

    _CACHE[key] = ""
    return ""
