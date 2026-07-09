"""ClinicalTrials.gov — clinical trial data (public domain, 商用利用可).

進行中・完了済み試験から競合状況・有効性シグナルを取得。
API v2: https://clinicaltrials.gov/api/v2/studies
"""
import requests

CT_API = "https://clinicaltrials.gov/api/v2/studies"

STATUS_LABEL = {
    "RECRUITING":              "Recruiting",
    "ACTIVE_NOT_RECRUITING":   "Active (not recruiting)",
    "COMPLETED":               "Completed",
    "TERMINATED":              "Terminated",
    "WITHDRAWN":               "Withdrawn",
    "NOT_YET_RECRUITING":      "Not yet recruiting",
    "ENROLLING_BY_INVITATION": "Enrolling by invitation",
    "UNKNOWN":                 "Unknown",
}


def _search(params: dict, max_results: int) -> list[dict]:
    """1回分の studies 検索 → 整形済みレコードのリストを返す（失敗時は空）。"""
    params = {
        **params,
        "pageSize": min(max_results * 2, 100),
        "format":   "json",
        "fields":   "NCTId,BriefTitle,OverallStatus,Phase,StartDate,"
                    "Condition,InterventionName,InterventionType,"
                    "LeadSponsorName,CollaboratorName",
    }
    try:
        r = requests.get(CT_API, params=params, timeout=20)
        r.raise_for_status()
        studies = r.json().get("studies", [])
    except Exception:
        return []

    results = []
    for s in studies:
        proto  = s.get("protocolSection", {})
        ident  = proto.get("identificationModule", {})
        status = proto.get("statusModule", {})
        design = proto.get("designModule", {})
        conds  = proto.get("conditionsModule", {})
        arms   = proto.get("armsInterventionsModule", {})
        sponsor_mod = proto.get("sponsorCollaboratorsModule", {})

        nct_id = ident.get("nctId", "")
        phase_list = design.get("phases", [])
        phase = "/".join(phase_list) if phase_list else "N/A"

        interventions = [
            iv.get("name", "")
            for iv in (arms.get("interventions") or [])
            if iv.get("type", "") in ("DRUG", "BIOLOGICAL", "GENETIC", "")
        ][:5]

        lead_sponsor  = (sponsor_mod.get("leadSponsor") or {}).get("name", "")
        collaborators = [c.get("name", "") for c in (sponsor_mod.get("collaborators") or [])][:5]

        results.append({
            "nct_id":        nct_id,
            "title":         ident.get("briefTitle", ""),
            "status":        STATUS_LABEL.get(status.get("overallStatus", ""), status.get("overallStatus", "")),
            "phase":         phase,
            "start_date":    status.get("startDateStruct", {}).get("date", ""),
            "conditions":    (conds.get("conditions") or [])[:3],
            "interventions": interventions,
            "sponsor":       lead_sponsor,
            "collaborators": collaborators,
            "url":           f"https://clinicaltrials.gov/study/{nct_id}" if nct_id else "",
        })
    return results


def get_trials(
    gene_symbol: str,
    disease: str,
    max_results: int = 10,
    drug_names: list[str] | None = None,
) -> list[dict]:
    """Return clinical trials related to a gene-disease pair.

    Runs two searches and merges/dedupes the results, because a single
    field can't catch everything:
      1. query.term=gene_symbol — matches the gene against all searchable
         fields (title, outcome, eligibility, etc). Catches genetic/
         biomarker studies, but NOT drug trials, since a trial's
         intervention is almost always the drug's brand/code name
         (e.g. "verubecestat"), not the target gene symbol.
      2. query.intr=<drug_names OR'd> — if known drugs for this target are
         passed in (from ChEMBL/OpenTargets/DGIdb), search for them by name
         directly. For BACE1 x Alzheimer's disease this alone finds ~4x
         more real trials (verubecestat, elenbecestat, lanabecestat,
         atabecestat, ...) than the gene-symbol search does, because most
         BACE1 inhibitor trials never mention "BACE1" in any searchable
         field.
    Returns:
        [{nct_id, title, status, phase, start_date, conditions,
          interventions, url}]
    """
    by_nct: dict[str, dict] = {}

    for r in _search({"query.cond": disease, "query.term": gene_symbol}, max_results):
        by_nct[r["nct_id"]] = r

    drug_names = [d for d in (drug_names or []) if d and d.strip()]
    if drug_names:
        # Essie 構文の OR で一括検索（API呼び出し回数を抑える）。長すぎる場合は分割。
        CHUNK = 15
        for i in range(0, len(drug_names), CHUNK):
            chunk = drug_names[i:i + CHUNK]
            intr_query = " OR ".join(chunk)
            for r in _search({"query.cond": disease, "query.intr": intr_query}, max_results):
                by_nct[r["nct_id"]] = r

    results = list(by_nct.values())
    # 直近の試験を優先（開始日降順、日付不明は末尾）
    results.sort(key=lambda r: r.get("start_date") or "", reverse=True)
    return results[:max_results]
