"""HPO (Human Phenotype Ontology) collector.

Uses HPO annotation flat files (cached in memory) — no dependency on
hpo.jax.org or Monarch APIs, both of which may be blocked.

Data sources:
  1. phenotype.hpoa  (disease → HP terms)
       https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa
  2. phenotype_to_genes.txt  (HP term → gene symbols)
       https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt
  3. OpenTargets GraphQL  (MONDO/EFO → OMIM cross-refs)
"""
from __future__ import annotations

import io
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import requests

_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json"})

# ── annotation file URLs ──────────────────────────────────────────────────
_HPOA_URL  = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"
_P2G_URL   = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt"
_OT_GQL    = "https://api.platform.opentargets.org/api/v4/graphql"

# ── in-memory cache ───────────────────────────────────────────────────────
_lock       = threading.Lock()
_loaded     = False
# disease_id (OMIM:xxx / ORPHA:xxx) → [(hpo_id, hpo_name, freq), ...]
_dis2pheno: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
# hpo_id → [gene_symbol, ...]
_hpo2genes: dict[str, list[str]] = defaultdict(list)


def _load_annotations(timeout: int = 30):
    """Download and parse HPO annotation files into memory (once)."""
    global _loaded
    with _lock:
        if _loaded:
            return
        try:
            # ── phenotype.hpoa → disease → HP term mapping ─────────────
            r = _SESSION.get(_HPOA_URL, timeout=timeout)
            r.raise_for_status()
            for line in io.StringIO(r.text):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                db_id   = parts[0]  # e.g. OMIM:201910
                hpo_id  = parts[3]  # e.g. HP:0001234
                freq    = parts[7] if len(parts) > 7 else ""
                # hpo_name is not in this file; fill from p2g later
                _dis2pheno[db_id].append((hpo_id, "", freq))

            # ── phenotype_to_genes.txt → HP term → gene symbols ────────
            r2 = _SESSION.get(_P2G_URL, timeout=timeout)
            r2.raise_for_status()
            hpo_names: dict[str, str] = {}
            for line in io.StringIO(r2.text):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                # format: hpo_id  hpo_name  ncbi_gene_id  gene_symbol  disease_id
                if len(parts) < 4:
                    continue
                hpo_id   = parts[0]
                hpo_name = parts[1]
                gene_sym = parts[3]
                if hpo_id not in hpo_names:
                    hpo_names[hpo_id] = hpo_name
                if gene_sym:
                    _hpo2genes[hpo_id].append(gene_sym)

            # backfill hpo_name in _dis2pheno
            for db_id, entries in _dis2pheno.items():
                _dis2pheno[db_id] = [
                    (hpo_id, hpo_names.get(hpo_id, hpo_id), freq)
                    for hpo_id, _, freq in entries
                ]
            _loaded = True
        except Exception as e:
            # leave _loaded=False so next call retries
            raise RuntimeError(f"HPO annotation files unavailable: {e}") from e


# ── MONDO → OMIM via OpenTargets ──────────────────────────────────────────

_OT_XREF_QUERY = """
query($efoId: String!) {
  disease(efoId: $efoId) {
    id
    name
    dbXRefs
  }
}
"""

def _get_omim_ids(mondo_id: str) -> list[str]:
    """Return OMIM IDs (e.g. ['OMIM:201910']) for a MONDO/EFO disease ID."""
    eid = mondo_id.replace("_", ":")  # MONDO_0008728 → MONDO:0008728
    try:
        resp = _SESSION.post(
            _OT_GQL,
            json={"query": _OT_XREF_QUERY, "variables": {"efoId": eid}},
            timeout=20,
        )
        resp.raise_for_status()
        xrefs = resp.json().get("data", {}).get("disease", {}).get("dbXRefs") or []
        return [x for x in xrefs if x.startswith("OMIM:")]
    except Exception:
        return []


# ── Public API ─────────────────────────────────────────────────────────────

def get_disease_phenotypes_from_cache(
    disease_id: str,
    max_phenotypes: int = 50,
) -> list[dict]:
    """Return list of {hpo_id, name, frequency} for a disease.

    disease_id: OMIM:xxx or ORPHA:xxx.
    """
    entries = _dis2pheno.get(disease_id, [])[:max_phenotypes]
    return [{"hpo_id": hpo_id, "name": name, "frequency": freq}
            for hpo_id, name, freq in entries]


def get_term_genes_from_cache(hpo_id: str) -> list[str]:
    """Return gene symbols associated with an HPO term."""
    return list(dict.fromkeys(_hpo2genes.get(hpo_id, [])))  # deduplicated


# ── Main evaluation function ───────────────────────────────────────────────

def evaluate_ppi_hpo_overlap(
    gene: str,
    ppi_partners: list[str],
    disease_name: str,
    omim_id: str = None,
    mondo_id: str = None,
    max_phenotypes: int = 30,
) -> dict:
    """Evaluate overlap between PPI partners and HPO symptom-associated genes.

    Steps:
      1. Load HPO annotation files (cached after first call).
      2. Resolve disease to OMIM ID (from mondo_id via OpenTargets, or omim_id directly).
      3. Look up HPO phenotype terms for the disease.
      4. For each term, get associated genes from cache.
      5. Compute overlap with ppi_partners.
    """
    ppi_set = {p.upper() for p in ppi_partners if p}

    # Step 1: load annotation files
    try:
        _load_annotations()
    except RuntimeError as e:
        return {
            "error": f"HPOアノテーションファイルにアクセスできません。"
                     f"インターネット接続またはプロキシ設定を確認してください。({e})",
            "disease_name": disease_name,
        }

    # Step 2: resolve OMIM ID
    resolved_omim_ids: list[str] = []
    disease_label = disease_name

    if omim_id:
        resolved_omim_ids = [omim_id]
    elif mondo_id:
        resolved_omim_ids = _get_omim_ids(mondo_id)

    if not resolved_omim_ids:
        return {
            "error": f"OMIM IDが取得できませんでした (MONDO={mondo_id})",
            "disease_name": disease_label,
        }

    # Step 3: collect phenotypes from all OMIM IDs for this disease
    phenos_seen: dict[str, dict] = {}
    for oid in resolved_omim_ids:
        for p in get_disease_phenotypes_from_cache(oid, max_phenotypes=max_phenotypes):
            phenos_seen[p["hpo_id"]] = p

    phenotypes = list(phenos_seen.values())[:max_phenotypes]

    if not phenotypes:
        return {
            "error": f"HPO症状データなし (OMIM={resolved_omim_ids})",
            "disease_name": disease_label,
            "disease_id": resolved_omim_ids[0] if resolved_omim_ids else "",
        }

    # Step 4: per-term overlap
    per_term = []
    all_hpo_genes: set[str] = set()
    for pheno in phenotypes:
        hpo_id = pheno["hpo_id"]
        tgenes_upper = {g.upper() for g in get_term_genes_from_cache(hpo_id) if g}
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

    # Step 5: aggregate
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

    return {
        "disease_id":     resolved_omim_ids[0] if resolved_omim_ids else "",
        "disease_name":   disease_label,
        "hpo_term_count": len(phenotypes),
        "per_term":       per_term,
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
