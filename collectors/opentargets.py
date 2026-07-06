"""OpenTargets Platform GraphQL API (Apache 2.0).

スキーマバージョン: v4（2024-2025）
主な変更点:
  knownDrugs → drugAndClinicalCandidates
  maximumClinicalTrialPhase → maximumClinicalStage
"""
import requests

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"

TARGET_QUERY = """
query($ensgId: String!) {
  target(ensemblId: $ensgId) {
    approvedSymbol
    approvedName
    biotype
    functionDescriptions
    associatedDiseases(enableIndirect: true, page: {index: 0, size: 20}) {
      rows {
        disease { id name }
        score
        datatypeScores { id score }
      }
    }
    drugAndClinicalCandidates {
      rows {
        drug { id name maximumClinicalStage }
        maxClinicalStage
        diseases { disease { name } }
      }
    }
  }
}
"""

SCORE_QUERY = """
query($efoId: String!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(enableIndirect: true, page: {index: 0, size: 500}) {
      rows {
        target { id approvedSymbol }
        score
        datatypeScores { id score }
      }
    }
  }
}
"""

SEARCH_QUERY = """
query($q: String!, $entity: [String!]) {
  search(queryString: $q, entityNames: $entity, page: {index: 0, size: 5}) {
    hits { id name entity }
  }
}
"""


def get_target_disease_evidence(
    gene_symbol: str,
    disease_name: str,
    gene_id: str = None,
    disease_id: str = None,
) -> dict:
    """Fetch OpenTargets evidence for a gene–disease pair."""

    # disease_id 解決
    disease_label = disease_name
    if not disease_id:
        r = requests.post(OT_API, json={
            "query": SEARCH_QUERY,
            "variables": {"q": disease_name, "entity": ["disease"]}
        }, timeout=20)
        r.raise_for_status()
        hits = [h for h in r.json().get("data", {}).get("search", {}).get("hits", [])
                if h.get("entity") == "disease"]
        if not hits:
            return {"error": f"Disease not found: {disease_name}"}
        disease_id    = hits[0]["id"]
        disease_label = hits[0]["name"]

    # gene_id 解決
    ensg_id = gene_id
    if not ensg_id:
        r = requests.post(OT_API, json={
            "query": SEARCH_QUERY,
            "variables": {"q": gene_symbol, "entity": ["target"]}
        }, timeout=20)
        r.raise_for_status()
        hits = [h for h in r.json().get("data", {}).get("search", {}).get("hits", [])
                if h.get("entity") == "target"]
        exact = [h for h in hits if h.get("name", "").upper() == gene_symbol.upper()]
        hits  = exact or hits
        ensg_id = hits[0]["id"] if hits else None

    # Target 情報 + 薬剤 + 関連疾患
    target_data = {}
    if ensg_id:
        r = requests.post(OT_API, json={
            "query": TARGET_QUERY,
            "variables": {"ensgId": ensg_id}
        }, timeout=25)
        r.raise_for_status()
        target_data = r.json().get("data", {}).get("target") or {}

    # 指定疾患に対するスコア
    assoc_score   = None
    datatype_scores = {}
    if ensg_id:
        try:
            r = requests.post(OT_API, json={
                "query": SCORE_QUERY,
                "variables": {"efoId": disease_id}
            }, timeout=30)
            r.raise_for_status()
            rows = (r.json().get("data", {}).get("disease") or {}) \
                       .get("associatedTargets", {}).get("rows", [])
            for row in rows:
                if (row.get("target") or {}).get("id") == ensg_id:
                    assoc_score     = row.get("score")
                    datatype_scores = {d["id"]: d["score"] for d in row.get("datatypeScores", [])}
                    break
        except Exception:
            pass

    # 薬剤リスト整形
    known_drugs = []
    dc = target_data.get("drugAndClinicalCandidates") or {}
    for row in (dc.get("rows") or []):
        if not isinstance(row, dict):
            continue
        drug = row.get("drug") or {}
        # "disease" フィールドが null の場合に備えて `or {}` でガード
        diseases = [
            (d.get("disease") or {}).get("name", "")
            for d in (row.get("diseases") or [])
            if isinstance(d, dict)
        ]
        known_drugs.append({
            "drug":      drug.get("name", ""),
            "max_phase": row.get("maxClinicalStage") or drug.get("maximumClinicalStage"),
            "disease":   diseases[0] if diseases else "",
            "mechanism": "",
        })

    # 関連疾患リスト整形
    assoc_dis_rows = (target_data.get("associatedDiseases") or {}).get("rows") or []
    associated_diseases = []
    for row in assoc_dis_rows:
        if not isinstance(row, dict):
            continue
        associated_diseases.append({
            "disease":         (row.get("disease") or {}).get("name", ""),
            "disease_id":      (row.get("disease") or {}).get("id", ""),
            "score":           row.get("score"),
            "datatype_scores": {d["id"]: d["score"] for d in (row.get("datatypeScores") or [])},
        })

    return {
        "gene_symbol":        gene_symbol,
        "ensembl_id":         ensg_id,
        "disease_id":         disease_id,
        "disease_label":      disease_label,
        "association_score":  assoc_score,
        "datatype_scores":    datatype_scores,
        "gene_info": {
            "name":     target_data.get("approvedName", ""),
            "biotype":  target_data.get("biotype", ""),
            "function": (target_data.get("functionDescriptions") or [""])[:2],
        },
        "known_drugs":         known_drugs,
        "associated_diseases": associated_diseases,
    }
