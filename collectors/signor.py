"""SIGNOR — Signaling Network Open Resource (CC BY 4.0, 商用利用可).

シグナル伝達ネットワークの因果関係データ（リン酸化・活性化・抑制など）。
APIキー不要。全データをTSVで取得し、ローカルにキャッシュして対象遺伝子でフィルタリング。
"""
import io
import time
import requests
from pathlib import Path

SIGNOR_TSV = "https://signor.uniroma2.it/getData.php?organism=9606&format=tsv"

_CACHE_DIR = Path(__file__).parent.parent / "ppi_cache"
_SIGNOR_CACHE = _CACHE_DIR / "signor_9606.tsv"
_CACHE_TTL = 7 * 24 * 3600  # 7日間


def _get_signor_tsv() -> str:
    """ローカルキャッシュから読み込む（古い場合は再ダウンロード）。"""
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if _SIGNOR_CACHE.exists():
        age = time.time() - _SIGNOR_CACHE.stat().st_mtime
        if age < _CACHE_TTL:
            return _SIGNOR_CACHE.read_text(encoding="utf-8")
    print("  [SIGNOR] TSV ダウンロード中（初回 or 7日経過）...")
    r = requests.get(SIGNOR_TSV, timeout=60)
    r.raise_for_status()
    _SIGNOR_CACHE.write_text(r.text, encoding="utf-8")
    return r.text

COLS = [
    "entityA", "typeA", "idA", "dbA",
    "entityB", "typeB", "idB", "dbB",
    "effect", "mechanism", "residue", "sequence",
    "taxId", "cellData", "tissueData", "modA", "modB",
    "pmid", "direct", "sentence_id", "annotated_by",
    "notes", "signor_id", "score",
]


def get_interactions(gene_symbol: str) -> list[dict]:
    """Return SIGNOR causal interactions involving the gene (as entityA or entityB)."""
    text = _get_signor_tsv()

    results = []
    gene_upper = gene_symbol.upper()

    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 9:
            continue

        entity_a = parts[0].strip().upper()
        entity_b = parts[4].strip().upper()

        if entity_a != gene_upper and entity_b != gene_upper:
            continue

        # タンパク質 or 複合体のみ
        if parts[1].strip() not in ("protein", "complex") and \
           parts[5].strip() not in ("protein", "complex"):
            continue

        partner    = parts[4].strip() if entity_a == gene_upper else parts[0].strip()
        direction  = "→" if entity_a == gene_upper else "←"
        effect     = parts[8].strip()
        mechanism  = parts[9].strip() if len(parts) > 9 else ""
        residue    = parts[10].strip() if len(parts) > 10 else ""
        pmid       = parts[17].strip() if len(parts) > 17 else ""
        score      = parts[23].strip() if len(parts) > 23 else ""

        try:
            score_f = float(score) if score else None
        except ValueError:
            score_f = None

        results.append({
            "source":    gene_symbol if entity_a == gene_upper else partner,
            "target":    partner if entity_a == gene_upper else gene_symbol,
            "partner":   partner,
            "direction": direction,
            "effect":    effect,
            "mechanism": mechanism,
            "residue":   residue,
            "pmid":      pmid,
            "score":     score_f,
            "db":        "SIGNOR",
        })

    return results
