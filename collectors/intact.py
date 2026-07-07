"""IntAct protein interaction database (EMBL-EBI, CC BY 4.0).

String代替として使用。商用利用可。
遺伝子ごとの結果を ppi_cache/intact/ にJSONキャッシュ（3日間有効）。

API 仕様変更 (2025):
  - participants フィールド削除 → moleculeA/B で相互作用分子を返す
  - intactScore → intactMiscore
  - publications list → publicationPubmedIdentifier (string)
  - detectionMethod/interactionType は文字列として返される
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
    center = gene_symbol.upper()
    interactions = []

    for item in data.get("content", []):
        mol_a = (item.get("moleculeA") or "").strip()
        mol_b = (item.get("moleculeB") or "").strip()
        type_a = (item.get("typeA") or "").strip().lower()
        type_b = (item.get("typeB") or "").strip().lower()

        # 相手分子（クエリ遺伝子でない方）と、その型を取得
        if mol_a.upper() == center:
            partner, partner_type = mol_b, type_b
        elif mol_b.upper() == center:
            partner, partner_type = mol_a, type_a
        else:
            # intactName でフォールバック
            name_a = (item.get("intactNameA") or "").upper()
            if center in name_a:
                partner, partner_type = mol_b, type_b
            else:
                partner, partner_type = mol_a, type_a

        if not partner:
            continue

        pubmed_id = item.get("publicationPubmedIdentifier", "")

        interactions.append({
            "interaction_id":  item.get("ac", ""),
            "partners":        [partner],
            # protein/peptide は遺伝子、small molecule 等はそれ以外として区別
            "partner_type":    "gene" if partner_type in ("protein", "peptide") else (partner_type or "unknown"),
            "detection_method": item.get("detectionMethod", ""),
            "interaction_type": item.get("type", ""),
            "confidence":      item.get("intactMiscore"),
            "pubmed_ids":      [pubmed_id] if pubmed_id else [],
        })

    _save_cache(gene_symbol, interactions)
    return interactions


def get_top_interactors(gene_symbol: str, top_n: int = 10) -> list[str]:
    """Return list of top interactor gene names."""
    interactions = get_interactions(gene_symbol, max_results=50)
    partners = []
    for ix in interactions:
        partners.extend(ix["partners"])
    from collections import Counter
    counts = Counter(partners)
    return [name for name, _ in counts.most_common(top_n)]
