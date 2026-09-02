"""HPO (Human Phenotype Ontology) collector.

Strategy (disease_gene_network01 method):
  1. OpenTargets `phenotypes` field: MONDO → HPO terms directly (primary)
  2. OpenTargets `dbXRefs` + HPO annotation file (fallback when OT phenotypes empty)
  3. HPO JAX API by disease ID (fallback)
  4. HPO JAX API search by name (fallback)

Genes per HP term:
  1. OpenTargets `associatedTargets` on HP term as efoId (primary)
  2. phenotype_to_genes.txt bulk file (fallback)

Scoring (disease_gene_network01 method):
  - Per-symptom hypergeometric test: P(X >= k) with M=20000, n=symptom_genes, N=ppi_count+1
  - Fisher's method: chi2 = -2 * sum(ln(p_i)), combined p from chi2 distribution
"""
from __future__ import annotations

import io
import math
import threading
from collections import defaultdict

import requests

_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json"})

_OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"
_HPO_BASE = "https://hpo.jax.org/api/hpo"

# ── OT: MONDO → HPO phenotypes (primary source) ──────────────────────────

_OT_PHENOTYPES_Q = """
query($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    phenotypes(page: { index: 0, size: $size }) {
      count
      rows {
        phenotypeHPO { id name description }
        evidence { aspect frequency qualifierNot resource }
      }
    }
  }
}
"""

_OT_HP_GENES_Q = """
query($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    associatedTargets(page: { index: 0, size: $size }) {
      rows {
        target { approvedSymbol }
        score
      }
    }
  }
}
"""

_OT_XREF_Q = """
query($efoId: String!) {
  disease(efoId: $efoId) {
    id
    name
    dbXRefs
  }
}
"""


def _ot_post(query: str, variables: dict, timeout: int = 25) -> dict:
    resp = _SESSION.post(_OT_GQL, json={"query": query, "variables": variables}, timeout=timeout)
    resp.raise_for_status()
    data = resp.json()
    if "errors" in data:
        raise RuntimeError(f"OT GraphQL: {data['errors'][0].get('message','')}")
    return data.get("data", {})


def _ot_get_disease_phenotypes(mondo_id: str, size: int = 50) -> tuple[str, list[dict]]:
    """Return (disease_name, [{hpo_id, name, frequency}, ...]) via OT phenotypes field.

    OT accepts MONDO_XXXXXXX or EFO_XXXXXXX directly (underscore format).
    Returns only non-excluded phenotypic abnormality terms (aspect P).
    """
    data = _ot_post(_OT_PHENOTYPES_Q, {"efoId": mondo_id, "size": size})
    dis = data.get("disease") or {}
    name = dis.get("name", "")
    pheno_data = dis.get("phenotypes") or {}
    rows = pheno_data.get("rows") or []
    result = []
    for row in rows:
        hpo_node = row.get("phenotypeHPO") or {}
        hpo_id = hpo_node.get("id", "")   # HP:0001234
        hpo_name = hpo_node.get("name", "")
        if not hpo_id or not hpo_name:
            continue
        evidences = row.get("evidence") or []
        # Skip if all evidence rows are excluded (qualifierNot=true)
        if evidences and all(e.get("qualifierNot") for e in evidences):
            continue
        # Extract frequency from first non-excluded evidence with a frequency value
        freq = ""
        for e in evidences:
            if not e.get("qualifierNot") and e.get("frequency"):
                freq = e["frequency"]
                break
        result.append({"hpo_id": hpo_id, "name": hpo_name, "frequency": freq})
    return name, result


def _ot_get_hpo_term_genes(hp_id: str, size: int = 200) -> list[str]:
    """Get gene symbols associated with an HP term via OT associatedTargets.

    OT uses underscore format: HP:0001234 → HP_0001234
    """
    ot_id = hp_id.replace(":", "_")  # HP:0001234 → HP_0001234
    try:
        data = _ot_post(_OT_HP_GENES_Q, {"efoId": ot_id, "size": size}, timeout=20)
        dis = data.get("disease") or {}
        rows = (dis.get("associatedTargets") or {}).get("rows") or []
        return [r["target"]["approvedSymbol"] for r in rows if r.get("target", {}).get("approvedSymbol")]
    except Exception:
        return []


def _ot_get_disease_xrefs(mondo_id: str) -> tuple[str, list[str]]:
    """Return (disease_name, [OMIM:xxx, ORPHA:xxx, ...]) via OT dbXRefs."""
    data = _ot_post(_OT_XREF_Q, {"efoId": mondo_id})
    dis = data.get("disease") or {}
    xrefs = dis.get("dbXRefs") or []
    name = dis.get("name", "")
    omim_ids = [x for x in xrefs if x.startswith("OMIM:") or x.startswith("ORPHA:")]
    return name, omim_ids


# ── HPO JAX API (fallback) ────────────────────────────────────────────────

def _hpo_parse_disease_data(data: dict) -> list[dict]:
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


def _hpo_api_disease_by_id(disease_id: str) -> list[dict]:
    """HPO JAX API: colon format (OMIM:201910, MONDO:0008728)."""
    uid = disease_id.replace("_", ":")
    r = _SESSION.get(f"{_HPO_BASE}/disease/{uid}", timeout=20)
    r.raise_for_status()
    return _hpo_parse_disease_data(r.json())


def _hpo_api_search_disease(query: str) -> list[dict]:
    r = _SESSION.get(f"{_HPO_BASE}/search",
                     params={"q": query, "max": 3, "category": "diseases"},
                     timeout=15)
    r.raise_for_status()
    diseases = r.json().get("diseases") or []
    if not diseases:
        return []
    disease_id = diseases[0].get("diseaseId", "")
    if not disease_id:
        return []
    r2 = _SESSION.get(f"{_HPO_BASE}/disease/{disease_id}", timeout=20)
    r2.raise_for_status()
    return _hpo_parse_disease_data(r2.json())


# ── HPO annotation files (last fallback) ─────────────────────────────────

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
_dis2pheno: dict[str, list] = defaultdict(list)
_hpo2genes: dict[str, list] = defaultdict(list)
_load_errors: list[str] = []


def _fetch_first(urls: list[str], timeout: int = 40) -> str | None:
    for url in urls:
        try:
            r = _SESSION.get(url, timeout=timeout)
            r.raise_for_status()
            return r.text
        except Exception as e:
            _load_errors.append(f"{url}: {e}")
    return None


def _ensure_hpoa():
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
            db_id  = parts[0]   # OMIM:201910 or ORPHA:xxx
            hpo_id = parts[3]   # HP:0001234
            freq   = parts[7] if len(parts) > 7 else ""
            _dis2pheno[db_id].append((hpo_id, "", freq))


def _ensure_p2g():
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
            if len(parts) < 4:
                continue
            hpo_id, hpo_name, _, gene_sym = parts[0], parts[1], parts[2], parts[3]
            if hpo_id not in hpo_names:
                hpo_names[hpo_id] = hpo_name
            if gene_sym:
                _hpo2genes[hpo_id].append(gene_sym)
        for db_id, entries in _dis2pheno.items():
            _dis2pheno[db_id] = [
                (hpo_id, hpo_names.get(hpo_id, hpo_id), freq)
                for hpo_id, _, freq in entries
            ]


# ── Statistical helpers ───────────────────────────────────────────────────

def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_sf(k: int, M: int, n: int, N: int) -> float:
    """P(X >= k) under Hypergeometric(M, n, N) in log space."""
    if k <= 0:
        return 1.0
    log_denom = _log_comb(M, N)
    total = 0.0
    for x in range(k, min(n, N) + 1):
        lp = _log_comb(n, x) + _log_comb(M - n, N - x) - log_denom
        total += math.exp(lp)
    return min(1.0, max(0.0, total))


def _fisher_combine(p_values: list[float]) -> float:
    """Fisher's method: chi2 = -2*sum(ln(p_i)), df=2m."""
    if not p_values:
        return 1.0
    floored = [max(p, 1e-300) for p in p_values]
    chi2 = -2.0 * sum(math.log(p) for p in floored)
    m = len(p_values)
    x = chi2 / 2.0
    term = 1.0
    s = 1.0
    for i in range(1, m):
        term *= x / i
        s += term
    return min(1.0, max(0.0, math.exp(-x) * s))


# ── Main evaluation ───────────────────────────────────────────────────────

def evaluate_ppi_hpo_overlap(
    gene: str,
    ppi_partners: list[str],
    disease_name: str,
    omim_id: str = None,
    mondo_id: str = None,
    max_phenotypes: int = 30,
) -> dict:
    """Evaluate overlap between target gene's PPI partners and HPO symptom genes.

    Phenotype retrieval priority:
      1. OT phenotypes field (MONDO → HP terms directly)
      2. OT dbXRefs → HPO annotation file
      3. HPO JAX API by ID
      4. HPO JAX API search by name

    Gene-per-HP retrieval priority:
      1. OT associatedTargets on HP term
      2. phenotype_to_genes.txt bulk file
    """
    ppi_set = {p.upper() for p in ppi_partners if p}
    disease_label = disease_name
    steps: list[str] = []
    omim_ids: list[str] = []
    phenotypes: list[dict] = []

    # ── Step 1: OT phenotypes field (primary — works for any MONDO/EFO ID) ──
    if mondo_id:
        try:
            name, phenotypes = _ot_get_disease_phenotypes(mondo_id, size=max_phenotypes * 2)
            if name:
                disease_label = name
            if phenotypes:
                steps.append(f"OT phenotypes: {mondo_id} → {len(phenotypes)} 症状")
            else:
                steps.append(f"OT phenotypes: {mondo_id} → 0件 (OTにHP対応なし)")
        except Exception as e:
            steps.append(f"OT phenotypes失敗: {str(e)[:80]}")

    # ── Step 2: OT dbXRefs → annotation file (when OT phenotypes empty) ─────
    if not phenotypes and mondo_id:
        try:
            name, omim_ids = _ot_get_disease_xrefs(mondo_id)
            if name and not disease_label:
                disease_label = name
            steps.append(f"OT dbXRefs: {mondo_id} → {omim_ids}")
        except Exception as e:
            steps.append(f"OT dbXRefs失敗: {str(e)[:60]}")

        if omim_ids:
            _ensure_hpoa()
            _ensure_p2g()
            steps.append(f"phenotype.hpoa: {len(_dis2pheno)} diseases")
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

    # ── Step 3: HPO JAX API by ID (OMIM, then MONDO) ────────────────────────
    if not phenotypes:
        if omim_id:
            omim_ids = [omim_id]
        for did in (omim_ids or []) + ([mondo_id] if mondo_id else []):
            try:
                phenotypes = _hpo_api_disease_by_id(did)
                if phenotypes:
                    steps.append(f"HPO API (ID): {did} → {len(phenotypes)} 症状")
                    break
            except Exception as e:
                steps.append(f"HPO API (ID) 失敗 ({did}): {str(e)[:60]}")

    # ── Step 4: HPO JAX API search by name ──────────────────────────────────
    if not phenotypes and disease_name:
        short_name = " ".join(disease_name.split()[:5])
        for q in ([disease_name] if disease_name != short_name else []) + [short_name]:
            try:
                phenotypes = _hpo_api_search_disease(q)
                if phenotypes:
                    steps.append(f"HPO API (検索): '{q}' → {len(phenotypes)} 症状")
                    break
            except Exception as e:
                steps.append(f"HPO API (検索) 失敗 '{q}': {str(e)[:60]}")

    if not phenotypes:
        return {
            "error": f"HPO症状データを取得できませんでした ({'; '.join(steps)})",
            "disease_name": disease_label,
            "steps": steps,
        }

    phenotypes = phenotypes[:max_phenotypes]

    # ── Step 5: per-term gene lookup + hypergeometric scoring ────────────────
    # OT associatedTargets is tried first for each HP term; p2g file is fallback.
    _ensure_p2g()
    p2g_ok = bool(_hpo2genes)
    steps.append(f"phenotype_to_genes fallback: {'OK' if p2g_ok else 'NG'} ({len(_hpo2genes)} terms)")

    M = 20_000          # background genome size
    N = len(ppi_set) + 1  # PPI partners + target itself

    per_term = []
    all_hpo_genes: set[str] = set()
    p_values_for_fisher: list[float] = []
    ot_gene_hits = 0

    for item in phenotypes:
        hpo_id   = item.get("hpo_id", "")
        hpo_name = item.get("name", hpo_id)
        freq     = item.get("frequency", "")

        # Gene lookup: OT first, then p2g file
        genes_raw = _ot_get_hpo_term_genes(hpo_id)
        if genes_raw:
            ot_gene_hits += 1
        else:
            genes_raw = _hpo2genes.get(hpo_id, [])

        tgenes_upper = {g.upper() for g in genes_raw if g}
        all_hpo_genes.update(tgenes_upper)
        overlap = sorted(ppi_set & tgenes_upper)
        k = len(overlap)
        n = len(tgenes_upper)

        p_val = _hypergeom_sf(k, M, n, N) if (k > 0 and n > 0) else 1.0
        if k > 0:
            p_values_for_fisher.append(p_val)

        per_term.append({
            "hpo_id":         hpo_id,
            "name":           hpo_name or hpo_id,
            "frequency":      freq,
            "hpo_gene_count": n,
            "overlap_genes":  overlap,
            "overlap_count":  k,
            "p_value":        p_val,
        })

    steps.append(f"HP遺伝子: OT {ot_gene_hits}/{len(phenotypes)} 件, p2g fallback {len(phenotypes)-ot_gene_hits} 件")
    per_term.sort(key=lambda x: x["p_value"])

    symptom_p_value = _fisher_combine(p_values_for_fisher) if p_values_for_fisher else 1.0

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
    note = ""
    if not p2g_ok and ot_gene_hits == 0:
        note = (f"HPO症状 {len(phenotypes)} 件を取得しましたが、"
                "遺伝子データを取得できず重複計算はスキップしました。")

    disease_id_out = omim_ids[0] if omim_ids else (mondo_id or "")

    return {
        "disease_id":     disease_id_out,
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
            "symptom_p_value":   symptom_p_value,
            "target_in_hpo":     gene.upper() in all_hpo_genes,
            "top_genes":         top_genes,
        },
    }
