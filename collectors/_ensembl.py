"""Shared Ensembl gene symbol → ID resolver with simple in-process cache."""
import requests

_CACHE: dict[str, str] = {}
ENSEMBL_REST = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"


def resolve_ensg(gene_symbol: str) -> str:
    """Return Ensembl ENSG ID for a human gene symbol, or '' on failure."""
    key = gene_symbol.upper()
    if key in _CACHE:
        return _CACHE[key]
    try:
        r = requests.get(
            f"{ENSEMBL_REST}/{gene_symbol}",
            headers={"Content-Type": "application/json"},
            timeout=10,
        )
        r.raise_for_status()
        ensg = r.json().get("id", "")
        _CACHE[key] = ensg
        return ensg
    except Exception:
        return ""
