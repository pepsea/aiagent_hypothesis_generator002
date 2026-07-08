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


def get_trials(gene_symbol: str, disease: str, max_results: int = 10) -> list[dict]:
    """Return clinical trials related to a gene-disease pair.

    Searches by disease condition; filters for interventions mentioning the gene.
    Returns:
        [{nct_id, title, status, phase, start_date, conditions,
          interventions, url}]
    """
    params = {
        "query.cond": disease,
        "query.intr": gene_symbol,
        "pageSize":   min(max_results * 2, 20),
        "format":     "json",
        "fields":     "NCTId,BriefTitle,OverallStatus,Phase,StartDate,"
                      "Condition,InterventionName,InterventionType,"
                      "LeadSponsorName,CollaboratorName",
    }

    try:
        r = requests.get(CT_API, params=params, timeout=20)
        r.raise_for_status()
        studies = r.json().get("studies", [])
    except Exception as e:
        return []

    results = []
    for s in studies[:max_results]:
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
