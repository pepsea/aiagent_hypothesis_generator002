"""SIGNOR — Signaling Network Open Resource (CC BY 4.0, 商用利用可).

シグナル伝達ネットワークの因果関係データ（リン酸化・活性化・抑制など）。
APIキー不要。全データをTSVで取得し、対象遺伝子でフィルタリング。
"""
import io
import requests

SIGNOR_TSV = "https://signor.uniroma2.it/getData.php?organism=9606&format=tsv"

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
    r = requests.get(SIGNOR_TSV, timeout=30)
    r.raise_for_status()

    results = []
    gene_upper = gene_symbol.upper()

    for line in r.text.splitlines():
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
