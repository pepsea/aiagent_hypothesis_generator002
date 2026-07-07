"""STRING protein interaction database (CC BY 4.0, 商用利用可).

信頼度スコア付きの機能的・物理的相互作用。required_score で閾値指定可能。
API: https://string-db.org/api/json/network
遺伝子ごとの結果を ppi_cache/string/ にJSONキャッシュ（3日間有効）。
"""
import json
import time
import requests
from pathlib import Path

BASE = "https://string-db.org/api/json/network"

_CACHE_DIR = Path(__file__).parent.parent / "ppi_cache" / "string"
_CACHE_TTL = 3 * 24 * 3600  # 3日間


def _cache_path(gene_symbol: str, required_score: int) -> Path:
    return _CACHE_DIR / f"{gene_symbol.upper()}_{required_score}.json"


def _load_cache(gene_symbol: str, required_score: int):
    p = _cache_path(gene_symbol, required_score)
    if p.exists() and time.time() - p.stat().st_mtime < _CACHE_TTL:
        return json.loads(p.read_text(encoding="utf-8"))
    return None


def _save_cache(gene_symbol: str, required_score: int, data: list[dict]):
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    _cache_path(gene_symbol, required_score).write_text(
        json.dumps(data, ensure_ascii=False), encoding="utf-8")


def get_interactions(gene_symbol: str, required_score: int = 400,
                     species: int = 9606, max_results: int = 100) -> list[dict]:
    """Return STRING interactions for a gene.

    required_score: 0–1000 の信頼度閾値（400=中信頼, 700=高信頼）。
    STRING は全パートナーが遺伝子/タンパク質（化合物は含まない）。
    """
    cached = _load_cache(gene_symbol, required_score)
    if cached is not None:
        return cached

    try:
        r = requests.get(BASE, params={
            "identifiers":   gene_symbol,
            "species":       species,
            "required_score": required_score,
            "limit":         max_results,
        }, timeout=25)
        if r.status_code == 404:
            _save_cache(gene_symbol, required_score, [])
            return []
        r.raise_for_status()
        rows = r.json()
    except Exception:
        return []

    center = gene_symbol.upper()
    interactions = []
    for row in rows:
        a = (row.get("preferredName_A") or "").strip()
        b = (row.get("preferredName_B") or "").strip()
        partner = b if a.upper() == center else a
        if not partner or partner.upper() == center:
            continue
        interactions.append({
            "source":       gene_symbol,
            "target":       partner,
            "partner":      partner,
            "partner_type": "gene",   # STRING は遺伝子のみ
            "effect":       "",
            "mechanism":    "functional/physical association",
            "direction":    "—",
            # STRING score は 0-1（表示・重み付け用）
            "score":        row.get("score"),
            # サブスコア（実験/データベース/共起など）
            "subscores": {
                "experimental": row.get("escore"),
                "database":     row.get("dscore"),
                "textmining":   row.get("tscore"),
                "coexpression": row.get("ascore"),
                "neighborhood": row.get("nscore"),
                "fusion":       row.get("fscore"),
                "cooccurrence": row.get("pscore"),
            },
            "pmid":         "",
            "db":           "STRING",
        })

    interactions.sort(key=lambda x: x.get("score") or 0, reverse=True)
    _save_cache(gene_symbol, required_score, interactions)
    return interactions
