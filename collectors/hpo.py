"""HPO (Human Phenotype Ontology) collector.

Strategy (in order of reliability):
  1. OpenTargets GraphQL — disease.phenotypes → HP term IDs  (OT already works)
  2. HPO annotation files (phenotype_to_genes.txt) — HP term → genes
     Tried from multiple mirror URLs; cached in memory after first load.
  3. Graceful error if all sources fail.
"""
from __future__ import annotations

import io
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json"})

# ── OpenTargets ───────────────────────────────────────────────────────────
_OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

_OT_PHENOTYPES_Q = """
query($efoId: String!) {
  disease(efoId: $efoId) {
    id
    name
    dbXRefs
    phenotypes {
      rows {
        phenotypeHPO { id name }
        frequencyHPO { id label }
      }
    }
  }
}
"""


def _ot_get_phenotypes(mondo_id: str) -> tuple[str, list[dict]]:
    """Query OT for HPO phenotype terms of a disease.

    Returns (omim_id_or_mondo, [{hpo_id, name, frequency}, ...])
    """
    eid = mondo_id.replace("_", ":")  # MONDO_0018479 → MONDO:0018479
    resp = _SESSION.post(
        _OT_GQL,
        json={"query": _OT_PHENOTYPES_Q, "variables": {"efoId": eid}},
        timeout=25,
    )
    resp.raise_for_status()
    dis = resp.json().get("data", {}).get("disease") or {}
    rows = (dis.get("phenotypes") or {}).get("rows") or []
    phenotypes = []
    for row in rows:
        hp = row.get("phenotypeHPO") or {}
        freq = row.get("frequencyHPO") or {}
        hpo_id = hp.get("id", "")
        name   = hp.get("name", "")
        freq_l = freq.get("label", "")
        if hpo_id and name:
            phenotypes.append({"hpo_id": hpo_id, "name": name, "frequency": freq_l})
    # also grab OMIM xref for annotation file lookup
    xrefs = dis.get("dbXRefs") or []
    omim_ids = [x for x in xrefs if x.startswith("OMIM:")]
    resolved_id = omim_ids[0] if omim_ids else eid
    return resolved_id, phenotypes


# ── HPO annotation file (HP term → genes) ────────────────────────────────
# Multiple candidate URLs — first that responds wins.
_P2G_URLS = [
    # purl canonical
    "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt",
    # GitHub releases (latest)
    "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype_to_genes.txt",
    # raw main branch (older location)
    "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/phenotype_to_genes.txt",
]

_lock        = threading.Lock()
_p2g_loaded  = False
_p2g_error   = ""
# hpo_id → [gene_symbol, ...]
_hpo2genes: dict[str, list[str]] = defaultdict(list)
# disease OMIM/ORPHA → [(hpo_id, name, freq), ...]
_dis2pheno: dict[str, list[tuple]] = defaultdict(list)


def _try_load_p2g() -> bool:
    """Try each URL in _P2G_URLS and parse phenotype_to_genes.txt."""
    global _p2g_error
    for url in _P2G_URLS:
        try:
            r = _SESSION.get(url, timeout=30)
            r.raise_for_status()
            hpo_names: dict[str, str] = {}
            for line in io.StringIO(r.text):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                # columns: hpo_id  hpo_name  ncbi_gene_id  gene_symbol  disease_id
                if len(parts) < 4:
                    continue
                hpo_id, hpo_name, _, gene_sym = parts[0], parts[1], parts[2], parts[3]
                if hpo_id not in hpo_names:
                    hpo_names[hpo_id] = hpo_name
                if gene_sym:
                    _hpo2genes[hpo_id].append(gene_sym)
            _p2g_error = ""
            return True
        except Exception as e:
            _p2g_error = str(e)
            continue
    return False


def _try_load_hpoa() -> bool:
    """Try each URL for phenotype.hpoa to build disease→pheno map."""
    hpoa_urls = [
        "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa",
        "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa",
    ]
    hpo_names: dict[str, str] = {k: "" for k in _hpo2genes}  # may already be populated
    for url in hpoa_urls:
        try:
            r = _SESSION.get(url, timeout=30)
            r.raise_for_status()
            for line in io.StringIO(r.text):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                db_id  = parts[0]   # OMIM:201910 / ORPHA:xxx
                hpo_id = parts[3]   # HP:xxxxxxx
                freq   = parts[7] if len(parts) > 7 else ""
                _dis2pheno[db_id].append((hpo_id, "", freq))
            return True
        except Exception:
            continue
    return False


def ensure_p2g_loaded():
    """Load phenotype_to_genes.txt into _hpo2genes (once, thread-safe)."""
    global _p2g_loaded
    with _lock:
        if not _p2g_loaded:
            _try_load_p2g()
            _p2g_loaded = True  # mark done even if failed (avoid repeated downloads)


# ── Public API ─────────────────────────────────────────────────────────────

def evaluate_ppi_hpo_overlap(
    gene: str,
    ppi_partners: list[str],
    disease_name: str,
    omim_id: str = None,
    mondo_id: str = None,
    max_phenotypes: int = 30,
) -> dict:
    """Evaluate overlap between PPI partners and HPO symptom-associated genes.

    Phenotype source priority:
      1. OpenTargets disease.phenotypes (via mondo_id)
      2. HPO annotation file _dis2pheno (via omim_id)

    Gene-per-term source: HPO annotation file _hpo2genes.
    """
    ppi_set = {p.upper() for p in ppi_partners if p}
    disease_label = disease_name

    # ── Step 1: get HP terms for the disease ─────────────────────────────
    phenotypes: list[dict] = []
    resolved_id = mondo_id or omim_id or ""

    if mondo_id:
        try:
            resolved_id, phenotypes = _ot_get_phenotypes(mondo_id)
            phenotypes = phenotypes[:max_phenotypes]
        except Exception as e:
            phenotypes = []
            resolved_id = mondo_id

    # fallback: annotation file disease map
    if not phenotypes and (omim_id or resolved_id):
        ensure_p2g_loaded()
        _try_load_hpoa()
        for oid in ([omim_id] if omim_id else []) + [resolved_id]:
            entries = _dis2pheno.get(oid, [])
            if entries:
                phenotypes = [
                    {"hpo_id": hid, "name": name, "frequency": freq}
                    for hid, name, freq in entries[:max_phenotypes]
                ]
                break

    if not phenotypes:
        return {
            "error": (
                "HPO症状データを取得できませんでした。"
                "OpenTargets API または HPO アノテーションファイルへの"
                f"アクセスを確認してください。(disease={mondo_id or omim_id})"
            ),
            "disease_name": disease_label,
        }

    # ── Step 2: load HP term → gene map ──────────────────────────────────
    ensure_p2g_loaded()

    # ── Step 3: per-term overlap ──────────────────────────────────────────
    per_term = []
    all_hpo_genes: set[str] = set()
    for pheno in phenotypes:
        hpo_id = pheno["hpo_id"]
        # deduplicated gene list from annotation file
        tgenes_upper = {g.upper() for g in _hpo2genes.get(hpo_id, []) if g}
        all_hpo_genes.update(tgenes_upper)
        overlap = sorted(ppi_set & tgenes_upper)
        per_term.append({
            "hpo_id":         hpo_id,
            "name":           pheno["name"],
            "frequency":      pheno.get("frequency", ""),
            "hpo_gene_count": len(tgenes_upper),
            "overlap_genes":  overlap,
            "overlap_count":  len(overlap),
        })

    per_term.sort(key=lambda x: x["overlap_count"], reverse=True)

    # ── Step 4: aggregate ─────────────────────────────────────────────────
    overlap_genes = ppi_set & all_hpo_genes
    gene_term_count: dict[str, int] = {}
    for pt in per_term:
        for g in pt["overlap_genes"]:
            gene_term_count[g.upper()] = gene_term_count.get(g.upper(), 0) + 1

    top_genes = sorted(
        [{"symbol": g, "term_count": c} for g, c in gene_term_count.items()],
        key=lambda x: x["term_count"], reverse=True,
    )[:15]

    total_hpo_genes = len(all_hpo_genes)
    overlap_score = round(
        len(overlap_genes) / max(1, min(len(ppi_set), total_hpo_genes)), 3
    ) if (ppi_set and total_hpo_genes) else 0.0

    # Note if gene map was empty (annotation file not loaded)
    note = ""
    if total_hpo_genes == 0 and phenotypes:
        note = (f"HPO症状 {len(phenotypes)} 件を取得しましたが、症状→遺伝子マッピングファイル"
                f"(phenotype_to_genes.txt)にアクセスできなかったため重複計算できません。")

    return {
        "disease_id":     resolved_id,
        "disease_name":   disease_label,
        "hpo_term_count": len(phenotypes),
        "per_term":       per_term,
        "note":           note,
        "summary": {
            "total_hpo_genes":   total_hpo_genes,
            "ppi_partner_count": len(ppi_set),
            "overlap_genes":     sorted(list(overlap_genes)),
            "overlap_count":     len(overlap_genes),
            "overlap_score":     overlap_score,
            "target_in_hpo":     gene.upper() in all_hpo_genes,
            "top_genes":         top_genes,
        },
    }
