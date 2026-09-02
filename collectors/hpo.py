"""HPO (Human Phenotype Ontology) collector.

Strategy (mirrors disease_gene_network01 exactly):

Phenotype (symptom) retrieval:
  1. OT `phenotypes` field — MONDO/EFO id → HP terms (multiple query variants, richest first)
  2. OT `dbXRefs` → phenotype.hpoa bulk file lookup by OMIM/ORPHA id
  3. phenotype.hpoa name search — when dbXRefs are empty, match disease name in bulk file
  4. HPO JAX API (last resort)

Genes per HP term (hybrid fetcher):
  1. OT `associatedTargets` on HP term as efoId (HP_XXXXXXX format)
  2. phenotype_to_genes.txt bulk file (fallback)

Scoring:
  - Per-symptom hypergeometric test: P(X >= k), M=20000, n=symptom_genes, N=ppi+1
  - Fisher's method combining all per-symptom p-values
"""
from __future__ import annotations

import io
import math
import re
import threading
from collections import defaultdict

import requests

_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json"})

_OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"
_HPO_BASE = "https://hpo.jax.org/api/hpo"

# ── OT phenotype queries — richest first, simpler fallback ───────────────
# Mirrors disease_gene_network01/collectors/opentargets.py PHENOTYPE_QUERIES
_OT_PHENOTYPE_QUERIES = [
    # Full: with evidence (frequency, qualifierNot, aspect, resource)
    """
query DiseasePhenotypes($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id name
    phenotypes(page: { index: 0, size: $size }) {
      count
      rows {
        phenotypeHPO { id name description }
        evidence { aspect frequency qualifierNot resource }
      }
    }
  }
}
""",
    # Without evidence details
    """
query DiseasePhenotypes($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id name
    phenotypes(page: { index: 0, size: $size }) {
      count
      rows { phenotypeHPO { id name description } }
    }
  }
}
""",
    # Minimal
    """
query DiseasePhenotypes($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id name
    phenotypes(page: { index: 0, size: $size }) {
      count
      rows { phenotypeHPO { id name } }
    }
  }
}
""",
]

_OT_HP_GENES_Q = """
query DiseaseTopGenes($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    associatedTargets(page: { index: 0, size: $size }) {
      rows { target { approvedSymbol } score }
    }
  }
}
"""

_OT_XREF_QUERIES = [
    "query DiseaseXrefs($efoId: String!) { disease(efoId: $efoId) { id name dbXRefs } }",
    "query DiseaseXrefs($efoId: String!) { disease(efoId: $efoId) { id name } }",
]

# Track which OT query variant works to avoid re-probing
_phenotype_query_index: int | None = None
_xref_query_index: int | None = None


def _ot_post(query: str, variables: dict, timeout: int = 25) -> dict:
    resp = _SESSION.post(_OT_GQL, json={"query": query, "variables": variables}, timeout=timeout)
    resp.raise_for_status()
    data = resp.json()
    if "errors" in data:
        raise RuntimeError(f"OT GraphQL: {data['errors'][0].get('message','unknown')}")
    return data.get("data", {})


def _ot_get_disease_phenotypes(mondo_id: str, size: int = 50) -> tuple[str, list[dict]]:
    """Return (disease_name, [{hpo_id, name, frequency}, ...]) via OT phenotypes field.

    Tries multiple query variants (richest→simplest) as disease_gene_network01 does.
    """
    global _phenotype_query_index
    order = ([_phenotype_query_index] if _phenotype_query_index is not None
             else range(len(_OT_PHENOTYPE_QUERIES)))

    disease_data = None
    for idx in order:
        try:
            data = _ot_post(_OT_PHENOTYPE_QUERIES[idx], {"efoId": mondo_id, "size": size})
            disease_data = data.get("disease")
            if disease_data is None:
                return "", []
            _phenotype_query_index = idx
            break
        except Exception:
            continue

    if disease_data is None:
        return "", []

    name = disease_data.get("name", "")
    rows = (disease_data.get("phenotypes") or {}).get("rows") or []
    result = []
    for row in rows:
        hpo_node = row.get("phenotypeHPO") or {}
        hpo_id = hpo_node.get("id", "")
        hpo_name = hpo_node.get("name", "")
        if not hpo_id or not hpo_name:
            continue
        evidence = [e for e in (row.get("evidence") or []) if isinstance(e, dict)]
        # Skip if ALL evidence rows are excluded (qualifierNot=True)
        if evidence and all(e.get("qualifierNot") for e in evidence):
            continue
        freq = next((e.get("frequency") for e in evidence
                     if e.get("frequency") and not e.get("qualifierNot")), "")
        result.append({"hpo_id": hpo_id, "name": hpo_name, "frequency": freq or ""})
    return name, result


def _normalise_xref(x: str) -> str:
    """Normalise a dbXRef to the form phenotype.hpoa uses.

    OT may return MIM:187300 or OMIM:187300; hpoa uses OMIM:187300.
    Orphanet appears as ORPHANET:, ORPHA:, ORPHACODE: → ORPHA:
    """
    if re.match(r"^MIM:\d", x, re.I):
        return "OMIM:" + x.split(":", 1)[1]
    if re.match(r"^(ORPHANET|ORPHACODE):", x, re.I):
        return "ORPHA:" + x.split(":", 1)[1]
    return x


def _ot_get_disease_xrefs(mondo_id: str) -> tuple[str, list[str]]:
    """Return (disease_name, [OMIM:xxx, ORPHA:xxx, ...]) — tries two query variants."""
    global _xref_query_index
    order = ([_xref_query_index] if _xref_query_index is not None
             else range(len(_OT_XREF_QUERIES)))

    for idx in order:
        try:
            data = _ot_post(_OT_XREF_QUERIES[idx], {"efoId": mondo_id})
            dis = data.get("disease") or {}
            name = dis.get("name", "")
            xrefs = dis.get("dbXRefs") or []
            _xref_query_index = idx
            omim_orpha = []
            for x in xrefs:
                norm = _normalise_xref(x)
                if re.match(r"^(OMIM|ORPHA):", norm, re.I):
                    omim_orpha.append(norm.upper())
            return name, omim_orpha
        except Exception:
            continue
    return "", []


def _ot_get_hpo_term_genes(hp_id: str, size: int = 200) -> list[str]:
    """Get gene symbols for an HP term via OT associatedTargets.

    HP:0001234 → HP_0001234 (OT uses underscore format).
    """
    ot_id = hp_id.replace(":", "_").upper()
    try:
        data = _ot_post(_OT_HP_GENES_Q, {"efoId": ot_id, "size": size}, timeout=20)
        rows = (data.get("disease") or {}).get("associatedTargets", {}).get("rows") or []
        return [r["target"]["approvedSymbol"] for r in rows
                if r.get("target", {}).get("approvedSymbol")]
    except Exception:
        return []


# ── HPO bulk annotation files ─────────────────────────────────────────────

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
# OMIM:xxx / ORPHA:xxx (upper) → [(hpo_id, disease_name, qualifier, freq, aspect), ...]
_dis2pheno: dict[str, list] = defaultdict(list)
# disease_name (lower) → [disease_id, ...]  (for name-based lookup)
_name2dis: dict[str, list[str]] = defaultdict(list)
# HP:xxx (upper) → [gene_symbol, ...]
_hpo2genes: dict[str, list] = defaultdict(list)
_hpo2names: dict[str, str] = {}
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
    """Load phenotype.hpoa → _dis2pheno + _name2dis (once)."""
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
            db_id    = parts[0].strip().upper()   # OMIM:201910
            dis_name = parts[1].strip() if len(parts) > 1 else ""
            qualifier= parts[2].strip() if len(parts) > 2 else ""
            hpo_id   = parts[3].strip()            # HP:0001234
            freq     = parts[7].strip() if len(parts) > 7 else ""
            aspect   = parts[10].strip() if len(parts) > 10 else ""
            if not db_id or not hpo_id.upper().startswith("HP:"):
                continue
            _dis2pheno[db_id].append((hpo_id, dis_name, qualifier, freq, aspect))
            if dis_name:
                _name2dis[dis_name.lower()].append(db_id)


def _ensure_p2g():
    """Load phenotype_to_genes.txt → _hpo2genes + _hpo2names (once)."""
    global _p2g_loaded
    with _lock:
        if _p2g_loaded:
            return
        _p2g_loaded = True
        text = _fetch_first(_P2G_URLS)
        if not text:
            return
        seen: dict[str, set] = defaultdict(set)
        for line in io.StringIO(text):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            hpo_id   = parts[0].strip().upper()
            hpo_name = parts[1].strip()
            gene_sym = parts[3].strip()
            if not hpo_id.startswith("HP:") or not gene_sym:
                continue
            if hpo_id not in _hpo2names and hpo_name:
                _hpo2names[hpo_id] = hpo_name
            if gene_sym.upper() not in seen[hpo_id]:
                seen[hpo_id].add(gene_sym.upper())
                _hpo2genes[hpo_id].append(gene_sym)


def _hpoa_phenotypes_for_ids(disease_ids: list[str], max_phenos: int) -> list[dict]:
    """Look up phenotype.hpoa by OMIM/ORPHA ids and return symptom list."""
    phenos_seen: dict[str, tuple] = {}
    for oid in disease_ids:
        key = oid.strip().upper()
        for entry in _dis2pheno.get(key, []):
            hpo_id, dis_name, qualifier, freq, aspect = entry
            hid = hpo_id.upper()
            # Skip non-phenotypic (inheritance etc.) and excluded
            if aspect and aspect.upper() not in ("P", ""):
                continue
            if qualifier.upper() == "NOT":
                continue
            if hid not in phenos_seen:
                phenos_seen[hid] = entry
    return [
        {"hpo_id": hid, "name": _hpo2names.get(hid, hid), "frequency": fr}
        for hid, (_, _, _, fr, _) in list(phenos_seen.items())[:max_phenos]
    ]


_STOPWORDS_NAME = {"the", "of", "to", "a", "an", "due", "by", "in", "for",
                   "with", "type", "form", "and", "or"}


def _name_words(text: str) -> set[str]:
    """Split a disease name into significant words (skip stopwords and short tokens)."""
    return {w for w in re.split(r"[\s,;/()]+", text.lower())
            if len(w) > 2 and w not in _STOPWORDS_NAME}


def _hpoa_find_ids_by_name(disease_name: str) -> list[str]:
    """Find OMIM/ORPHA ids by name matching in phenotype.hpoa.

    Strategy (mirrors disease_gene_network01):
      1. Exact case-insensitive match
      2. Word-set match: wanted words all found in candidate name (any order)
         This handles OMIM's inverted naming like "Telangiectasia, hereditary hemorrhagic"
    """
    wanted = (disease_name or "").strip().lower()
    if not wanted:
        return []
    # Exact match
    hits = list(dict.fromkeys(_name2dis.get(wanted, [])))
    if hits:
        return hits[:5]
    # Word-set match — handles inverted/parenthetical OMIM names
    wanted_words = _name_words(wanted)
    if not wanted_words:
        return hits[:5]
    scored: list[tuple[float, str]] = []
    for name_key, ids in _name2dis.items():
        cand_words = _name_words(name_key)
        if not cand_words:
            continue
        common = wanted_words & cand_words
        ratio = len(common) / len(wanted_words)
        if ratio >= 0.8:
            for did in ids:
                scored.append((ratio, did))
    scored.sort(reverse=True)
    seen: set[str] = set()
    for _, did in scored:
        if did not in seen:
            seen.add(did)
            hits.append(did)
            if len(hits) >= 5:
                break
    return hits


# ── HPO JAX API (last resort) ─────────────────────────────────────────────

def _hpo_parse_disease_data(data: dict) -> list[dict]:
    raw = data.get("catTermsCombo") or data.get("phenotypes") or []
    result = []
    for p in raw:
        hpo_id = p.get("ontologyId") or (p.get("term") or {}).get("id")
        name   = p.get("name") or (p.get("term") or {}).get("name")
        freq   = p.get("frequency") or {}
        freq_l = freq.get("label", "") if isinstance(freq, dict) else str(freq or "")
        if hpo_id and name:
            result.append({"hpo_id": hpo_id, "name": name, "frequency": freq_l})
    return result


def _hpo_api_disease_by_id(disease_id: str) -> list[dict]:
    uid = disease_id.replace("_", ":")
    r = _SESSION.get(f"{_HPO_BASE}/disease/{uid}", timeout=20)
    r.raise_for_status()
    return _hpo_parse_disease_data(r.json())


def _hpo_api_search_disease(query: str) -> list[dict]:
    r = _SESSION.get(f"{_HPO_BASE}/search",
                     params={"q": query, "max": 3, "category": "diseases"}, timeout=15)
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


# ── Statistical helpers ───────────────────────────────────────────────────

def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_sf(k: int, M: int, n: int, N: int) -> float:
    """P(X >= k) under Hypergeometric(M, n, N), log-space."""
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
    """Evaluate overlap between PPI partners and HPO symptom genes.

    Phenotype retrieval (disease_gene_network01 strategy):
      1. OT phenotypes field — direct MONDO→HP (multiple query variants)
      2. OT dbXRefs → phenotype.hpoa file
      3. phenotype.hpoa name search (when dbXRefs empty)
      4. HPO JAX API (last resort)

    Gene retrieval per HP term:
      1. OT associatedTargets on HP_XXXXXXX
      2. phenotype_to_genes.txt fallback
    """
    ppi_set = {p.upper() for p in ppi_partners if p}
    disease_label = disease_name
    steps: list[str] = []
    omim_ids: list[str] = []
    phenotypes: list[dict] = []

    # ── Step 1: OT phenotypes field ──────────────────────────────────────────
    if mondo_id:
        try:
            name, phenotypes = _ot_get_disease_phenotypes(mondo_id, size=max_phenotypes * 2)
            if name:
                disease_label = name
            steps.append(f"OT phenotypes: {mondo_id} → {len(phenotypes)} 症状")
        except Exception as e:
            steps.append(f"OT phenotypes失敗: {str(e)[:80]}")

    # ── Step 2: OT dbXRefs → phenotype.hpoa ─────────────────────────────────
    if not phenotypes:
        if omim_id:
            omim_ids = [omim_id]
            steps.append(f"OMIM直接指定: {omim_id}")
        elif mondo_id:
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
            phenotypes = _hpoa_phenotypes_for_ids(omim_ids, max_phenotypes)
            if phenotypes:
                steps.append(f"phenotype.hpoa (ID): {len(phenotypes)} 症状")
            else:
                steps.append(f"phenotype.hpoa (ID): 0件 ({len(_dis2pheno)} diseases indexed)")

    # ── Step 3: phenotype.hpoa name search (disease_gene_network01 fallback) ─
    if not phenotypes:
        _ensure_hpoa()
        _ensure_p2g()
        name_ids = _hpoa_find_ids_by_name(disease_label or disease_name)
        steps.append(f"phenotype.hpoa name search: '{disease_label}' → {name_ids}")
        if name_ids:
            phenotypes = _hpoa_phenotypes_for_ids(name_ids, max_phenotypes)
            if phenotypes:
                omim_ids = omim_ids or name_ids
                steps.append(f"phenotype.hpoa (name): {len(phenotypes)} 症状")

    # ── Step 4: HPO JAX API ──────────────────────────────────────────────────
    if not phenotypes:
        for did in (omim_ids or []) + ([mondo_id] if mondo_id else []):
            try:
                phenotypes = _hpo_api_disease_by_id(did)
                if phenotypes:
                    steps.append(f"HPO API (ID): {did} → {len(phenotypes)} 症状")
                    break
            except Exception as e:
                steps.append(f"HPO API (ID) 失敗 ({did}): {str(e)[:60]}")

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
    _ensure_p2g()
    p2g_ok = bool(_hpo2genes)

    M = 20_000
    N = len(ppi_set) + 1

    per_term = []
    all_hpo_genes: set[str] = set()
    p_values_for_fisher: list[float] = []
    ot_gene_hits = 0

    for item in phenotypes:
        hpo_id   = item.get("hpo_id", "")
        hpo_name = item.get("name", hpo_id)
        freq     = item.get("frequency", "")

        # Gene lookup: OT associatedTargets first, then p2g file
        genes_raw = _ot_get_hpo_term_genes(hpo_id)
        if genes_raw:
            ot_gene_hits += 1
        else:
            hpo_key = hpo_id.strip().upper()
            genes_raw = _hpo2genes.get(hpo_key, [])
            if not genes_raw and not hpo_key.startswith("HP:"):
                hpo_key = hpo_key.replace("_", ":")
                genes_raw = _hpo2genes.get(hpo_key, [])

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

    steps.append(f"HP遺伝子: OT {ot_gene_hits}/{len(phenotypes)}, p2g {len(phenotypes)-ot_gene_hits}")
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
        note = (f"HPO症状 {len(phenotypes)} 件取得しましたが、"
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
