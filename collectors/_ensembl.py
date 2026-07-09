"""Shared Ensembl gene symbol → ID resolver with simple in-process cache."""
import time

import requests

_CACHE: dict[str, str] = {}
ENSEMBL_REST = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"


def resolve_ensg(gene_symbol: str, max_retries: int = 3) -> str:
    """Return Ensembl ENSG ID for a human gene symbol, or '' on failure.

    404（遺伝子が存在しない）は即座に諦めるが、タイムアウト・5xx等の一時的な
    障害は短い待機を挟んでリトライする（HPA/GTEx など呼び出し元共通）。
    """
    key = gene_symbol.upper()
    if key in _CACHE:
        return _CACHE[key]

    last_exc = None
    for attempt in range(max_retries):
        try:
            r = requests.get(
                f"{ENSEMBL_REST}/{gene_symbol}",
                headers={"Content-Type": "application/json"},
                timeout=10,
            )
            if r.status_code == 404:
                return ""
            r.raise_for_status()
            ensg = r.json().get("id", "")
            _CACHE[key] = ensg
            return ensg
        except Exception as e:
            last_exc = e
            if attempt < max_retries - 1:
                time.sleep(1 + attempt * 2)
    return ""
