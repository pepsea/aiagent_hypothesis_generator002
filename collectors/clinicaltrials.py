"""ClinicalTrials.gov — clinical trial data (public domain, 商用利用可).

進行中・完了済み試験から競合状況・有効性シグナルを取得。
API v2: https://clinicaltrials.gov/api/v2/studies

試験件数そのものが競合状況の指標になるため（同じ標的・疾患を狙う治験が
多いほど新規参入のリスクが高い）、表示・取得件数に上限は設けず全件取得する。
"""
import requests

CT_API = "https://clinicaltrials.gov/api/v2/studies"
MAX_PAGES = 20  # 1ページ1000件 × 20 = 最大2万件（安全上限、通常は数ページで尽きる）

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


def _search(params: dict) -> list[dict]:
    """該当する studies を全ページ取得 → 整形済みレコードのリストを返す（失敗時は空）。"""
    base_params = {
        **params,
        "pageSize": 1000,
        "format":   "json",
        "fields":   "NCTId,BriefTitle,OverallStatus,Phase,StartDate,"
                    "Condition,InterventionName,InterventionType,"
                    "LeadSponsorName,CollaboratorName",
    }

    studies = []
    page_token = None
    try:
        for _ in range(MAX_PAGES):
            p = dict(base_params)
            if page_token:
                p["pageToken"] = page_token
            r = requests.get(CT_API, params=p, timeout=20)
            r.raise_for_status()
            data = r.json()
            studies.extend(data.get("studies", []))
            page_token = data.get("nextPageToken")
            if not page_token:
                break
    except Exception:
        pass  # ここまでに取得できた分は活かす

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
    drug_names: list[str] | None = None,
) -> list[dict]:
    """Return ALL clinical trials related to a gene-disease pair (no cap).

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

    件数の上限は設けない（試験数自体が競合の激しさ＝新規参入リスクの指標に
    なるため）。直近の試験を先頭にして返す。
    Returns:
        [{nct_id, title, status, phase, start_date, conditions,
          interventions, url}]
    """
    by_nct: dict[str, dict] = {}

    for r in _search({"query.cond": disease, "query.term": gene_symbol}):
        by_nct[r["nct_id"]] = r

    drug_names = [d for d in (drug_names or []) if d and d.strip()]
    if drug_names:
        # Essie 構文の OR で一括検索（API呼び出し回数を抑える）。長すぎる場合は分割。
        CHUNK = 15
        for i in range(0, len(drug_names), CHUNK):
            chunk = drug_names[i:i + CHUNK]
            intr_query = " OR ".join(chunk)
            for r in _search({"query.cond": disease, "query.intr": intr_query}):
                by_nct[r["nct_id"]] = r

    results = list(by_nct.values())
    # 直近の試験を優先（開始日降順、日付不明は末尾）
    results.sort(key=lambda r: r.get("start_date") or "", reverse=True)
    return results
