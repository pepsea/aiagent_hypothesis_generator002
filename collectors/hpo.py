"""HPO (Human Phenotype Ontology) collector.

Strategy:
  1. OpenTargets GraphQL: MONDO → OMIM cross-ref (dbXRefs)
  2. HPO annotation files: disease (OMIM) → phenotypes, phenotype → genes
     Tries purl.obolibrary.org and GitHub release URLs; cached in memory.
"""
from __future__ import annotations

import io
import threading
from collections import defaultdict

import requests

_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json"})

# ── OpenTargets: MONDO → OMIM/dbXRefs ───────────────────────────────────
_OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

_OT_XREF_Q = """
query($efoId: String!) {
  disease(efoId: $efoId) {
    id
    name
    dbXRefs
  }
}
"""

_HPO_BASE = "https://hpo.jax.org/api/hpo"


def _ot_get_disease_info(mondo_id: str) -> tuple[str, list[str]]:
    """Return (disease_name, [OMIM:xxx, ...]) via OpenTargets GraphQL."""
    eid = mondo_id.replace("_", ":")
    resp = _SESSION.post(
        _OT_GQL,
        json={"query": _OT_XREF_Q, "variables": {"efoId": eid}},
        timeout=25,
    )
    resp.raise_for_status()
    data = resp.json()
    if "errors" in data:
        raise RuntimeError(f"OT GraphQL error: {data['errors'][0].get('message','')}")
    dis = data.get("data", {}).get("disease") or {}
    xrefs = dis.get("dbXRefs") or []
    name  = dis.get("name", "")
    omim_ids = [x for x in xrefs if x.startswith("OMIM:")]
    return name, omim_ids


def _hpo_api_disease(disease_id: str) -> list[dict]:
    """Try HPO JAX API directly with any disease ID (OMIM:, ORPHA:, MONDO:).

    Returns [{hpo_id, name, frequency}, ...] or [] on failure.
    """
    uid = disease_id.replace("_", ":")  # MONDO_0008728 → MONDO:0008728
    r = _SESSION.get(f"{_HPO_BASE}/disease/{uid}", timeout=20)
    r.raise_for_status()
    data = r.json()
    raw = data.get("catTermsCombo") or data.get("phenotypes") or []
    result = []
    for p in raw:
        hpo_id = p.get("ontologyId") or (p.get("term") or {}).get("id")
        name   = p.get("name")        or (p.get("term") or {}).get("name")
        freq   = p.get("frequency") or {}
        freq_l = freq.get("label", "") if isinstance(freq, dict) else str(freq or "")
        if hpo_id and name:
            result.append({"hpo_id": hpo_id, "name": name, "frequency": freq_l})
    return result


# ── HPO annotation files ──────────────────────────────────────────────────
_HPOA_URLS = [
    "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa",
    "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa",
]
_P2G_URLS = [
    "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt",
    "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype_to_genes.txt",
]

_lock = threading.Lock()
_hpoa_loaded = False
_p2g_loaded  = False
# OMIM:xxx / ORPHA:xxx → [(hpo_id, hpo_name, freq_str), ...]
_dis2pheno: dict[str, list] = defaultdict(list)
# HP:xxx → [gene_symbol, ...]
_hpo2genes: dict[str, list] = defaultdict(list)
_load_errors: list[str] = []


def _fetch_first(urls: list[str], timeout: int = 40) -> str | None:
    """Download the first URL that succeeds; return text content or None."""
    for url in urls:
        try:
            r = _SESSION.get(url, timeout=timeout)
            r.raise_for_status()
            return r.text
        except Exception as e:
            _load_errors.append(f"{url}: {e}")
    return None


def _ensure_hpoa():
    """Load phenotype.hpoa → _dis2pheno (once)."""
    global _hpoa_loaded
    with _lock:
        if _hpoa_loaded:
            return
        _hpoa_loaded = True
        text = _fetch_first(_HPOA_URLS)
        if not text:
            return
        for line in io.StringIO(text):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            db_id  = parts[0]  # OMIM:201910
            hpo_id = parts[3]  # HP:0001234
            freq   = parts[7] if len(parts) > 7 else ""
            _dis2pheno[db_id].append((hpo_id, "", freq))


def _ensure_p2g():
    """Load phenotype_to_genes.txt → _hpo2genes (once)."""
    global _p2g_loaded
    with _lock:
        if _p2g_loaded:
            return
        _p2g_loaded = True
        text = _fetch_first(_P2G_URLS)
        if not text:
            return
        hpo_names: dict[str, str] = {}
        for line in io.StringIO(text):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            # hpo_id  hpo_name  ncbi_gene_id  gene_symbol  disease_id
            if len(parts) < 4:
                continue
            hpo_id, hpo_name, _, gene_sym = parts[0], parts[1], parts[2], parts[3]
            if hpo_id not in hpo_names:
                hpo_names[hpo_id] = hpo_name
            if gene_sym:
                _hpo2genes[hpo_id].append(gene_sym)
        # backfill hpo_name into _dis2pheno
        for db_id, entries in _dis2pheno.items():
            _dis2pheno[db_id] = [
                (hpo_id, hpo_names.get(hpo_id, hpo_id), freq)
                for hpo_id, _, freq in entries
            ]


# ── Main evaluation ────────────────────────────────────────────────────────

def evaluate_ppi_hpo_overlap(
    gene: str,
    ppi_partners: list[str],
    disease_name: str,
    omim_id: str = None,
    mondo_id: str = None,
    max_phenotypes: int = 30,
) -> dict:
    """Evaluate overlap between target gene's PPI partners and HPO symptom genes."""
    ppi_set = {p.upper() for p in ppi_partners if p}
    disease_label = disease_name
    steps: list[str] = []

    # Step 1a: OMIM IDs — direct or via OT dbXRefs
    omim_ids: list[str] = []
    if omim_id:
        omim_ids = [omim_id]
        steps.append(f"OMIM直接指定: {omim_id}")
    elif mondo_id:
        try:
            name, omim_ids = _ot_get_disease_info(mondo_id)
            if name:
                disease_label = name
            steps.append(f"OT dbXRefs: {mondo_id} → {omim_ids}")
        except Exception as e:
            steps.append(f"OT dbXRefs失敗: {e}")

    # Step 1b: try HPO API directly (MONDO or OMIM) — works if hpo.jax.org accessible
    phenotypes: list[dict] = []
    for did in (omim_ids or []) + ([mondo_id] if mondo_id else []):
        try:
            phenotypes = _hpo_api_disease(did)
            if phenotypes:
                steps.append(f"HPO API直接取得: {did} → {len(phenotypes)} 症状")
                break
        except Exception as e:
            steps.append(f"HPO API失敗 ({did}): {e}")

    # Step 1c: fallback — HPO annotation files
    if not phenotypes:
        _ensure_hpoa()
        _ensure_p2g()
        hpoa_ok = bool(_dis2pheno)
        p2g_ok  = bool(_hpo2genes)
        steps.append(f"phenotype.hpoa: {'OK' if hpoa_ok else 'NG'} ({len(_dis2pheno)} diseases)")

        phenos_seen: dict[str, tuple] = {}
        for oid in omim_ids:
            for entry in _dis2pheno.get(oid, []):
                hid = entry[0]
                if hid not in phenos_seen:
                    phenos_seen[hid] = entry
        if phenos_seen:
            phenotypes = [
                {"hpo_id": hid, "name": nm, "frequency": fr}
                for hid, nm, fr in list(phenos_seen.values())[:max_phenotypes]
            ]
            steps.append(f"annotation file: {len(phenotypes)} 症状")

    if not phenotypes:
        return {
            "error": f"HPO症状データを取得できませんでした ({'; '.join(steps)})",
            "disease_name": disease_label,
            "steps": steps,
        }

    phenotypes = phenotypes[:max_phenotypes]

    # Step 2: load HP→gene map if not already loaded
    if not _hpo2genes:
        _ensure_p2g()
    p2g_ok = bool(_hpo2genes)
    steps.append(f"phenotype_to_genes: {'OK' if p2g_ok else 'NG'} ({len(_hpo2genes)} terms)")

    # Step 3: per-term overlap
    per_term = []
    all_hpo_genes: set[str] = set()
    for hpo_id, hpo_name, freq in phenotypes:
        tgenes_upper = {g.upper() for g in _hpo2genes.get(hpo_id, []) if g}
        all_hpo_genes.update(tgenes_upper)
        overlap = sorted(ppi_set & tgenes_upper)
        per_term.append({
            "hpo_id":         hpo_id,
            "name":           hpo_name or hpo_id,
            "frequency":      freq,
            "hpo_gene_count": len(tgenes_upper),
            "overlap_genes":  overlap,
            "overlap_count":  len(overlap),
        })

    per_term.sort(key=lambda x: x["overlap_count"], reverse=True)

    # Step 4: aggregate
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

    note = ""
    if not p2g_ok:
        note = (f"HPO症状 {len(phenotypes)} 件を取得しましたが、"
                "phenotype_to_genes.txt にアクセスできず遺伝子重複計算はスキップしました。")

    return {
        "disease_id":     omim_ids[0],
        "disease_name":   disease_label,
        "hpo_term_count": len(phenotypes),
        "per_term":       per_term,
        "note":           note,
        "steps":          steps,
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
