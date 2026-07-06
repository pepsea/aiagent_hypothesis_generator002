"""IntAct protein interaction database (EMBL-EBI, CC BY 4.0).

String代替として使用。商用利用可。
遺伝子ごとの結果を ppi_cache/intact/ にJSONキャッシュ（3日間有効）。
"""
import json
import time
import requests
from pathlib import Path

BASE = "https://www.ebi.ac.uk/intact/ws/interaction"

_CACHE_DIR = Path(__file__).parent.parent / "ppi_cache" / "intact"
_CACHE_TTL = 3 * 24 * 3600  # 3日間


def _cache_path(gene_symbol: str) -> Path:
    return _CACHE_DIR / f"{gene_symbol.upper()}.json"


def _load_cache(gene_symbol: str) -> list[dict] | None:
    p = _cache_path(gene_symbol)
    if p.exists() and time.time() - p.stat().st_mtime < _CACHE_TTL:
        return json.loads(p.read_text(encoding="utf-8"))
    return None


def _save_cache(gene_symbol: str, data: list[dict]):
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    _cache_path(gene_symbol).write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")


def get_interactions(gene_symbol: str, species: int = 9606, max_results: int = 20) -> list[dict]:
    """Return top interactors for a gene from IntAct."""
    cached = _load_cache(gene_symbol)
    if cached is not None:
        return cached

    r = requests.get(f"{BASE}/findInteractions/{gene_symbol}", params={
        "page": 0, "pageSize": max_results,
        "query": f"species:{species}",
    }, timeout=20)

    if r.status_code == 404:
        _save_cache(gene_symbol, [])
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

    _save_cache(gene_symbol, interactions)
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
