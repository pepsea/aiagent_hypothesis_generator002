"""OpenTargets Platform GraphQL API (Apache 2.0).

スキーマバージョン: v4（2024-2025）
主な変更点:
  knownDrugs → drugAndClinicalCandidates
  maximumClinicalTrialPhase → maximumClinicalStage

enableIndirect: true でウェブサイト表示と一致させている。
オントロジー階層を辿った間接エビデンスを含むため、
platform.opentargets.org の表示スコアと同じ値が得られる。
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
    synonyms { terms scope }
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

# 疾患の上位関連遺伝子を取得（パスウェイ解析用）
DISEASE_TOP_GENES_QUERY = """
query($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    associatedTargets(enableIndirect: true, page: {index: 0, size: $size}) {
      rows {
        target { id approvedSymbol }
        score
      }
    }
  }
}
"""


def get_disease_top_genes(disease_id: str, top_n: int = 20) -> list[dict]:
    """疾患に関連するスコア上位の遺伝子リストを返す（パスウェイ解析用）。"""
    try:
        r = requests.post(OT_API, json={
            "query": DISEASE_TOP_GENES_QUERY,
            "variables": {"efoId": disease_id, "size": top_n},
        }, timeout=20)
        r.raise_for_status()
        rows = (r.json().get("data", {}).get("disease") or {}) \
                   .get("associatedTargets", {}).get("rows", [])
        return [
            {"symbol": row["target"]["approvedSymbol"], "score": row["score"]}
            for row in rows if row.get("target")
        ]
    except Exception:
        return []


# 遺伝子×疾患ペアを直接指定してスコアを取得（size制限の影響を受けない）
PAIR_SCORE_QUERY = """
query($ensgId: String!, $efoId: String!) {
  target(ensemblId: $ensgId) {
    approvedSymbol
    associatedDiseases(
      enableIndirect: true
      filter: { ids: [$efoId] }
      page: { index: 0, size: 1 }
    ) {
      rows {
        disease { id name }
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

    # 指定疾患に対するスコア（遺伝子×疾患ペアを直接指定して確実に取得）
    assoc_score      = None
    datatype_scores  = {}
    disease_synonyms = []
    if ensg_id and disease_id:
        # ① ペアクエリでスコアを取得（size制限の影響を受けない）
        try:
            r = requests.post(OT_API, json={
                "query": PAIR_SCORE_QUERY,
                "variables": {"ensgId": ensg_id, "efoId": disease_id}
            }, timeout=20)
            r.raise_for_status()
            pair_rows = (
                (r.json().get("data", {}).get("target") or {})
                .get("associatedDiseases", {})
                .get("rows", [])
            )
            if pair_rows:
                assoc_score     = pair_rows[0].get("score")
                datatype_scores = {
                    d["id"]: d["score"]
                    for d in pair_rows[0].get("datatypeScores", [])
                }
        except Exception:
            pass

        # ② synonyms は disease クエリから取得
        try:
            r = requests.post(OT_API, json={
                "query": SCORE_QUERY,
                "variables": {"efoId": disease_id}
            }, timeout=30)
            r.raise_for_status()
            disease_data = r.json().get("data", {}).get("disease") or {}
            for syn in disease_data.get("synonyms", []):
                disease_synonyms.extend(syn.get("terms", []))
            # ペアクエリでスコアが取れなかった場合のフォールバック
            if assoc_score is None:
                rows = disease_data.get("associatedTargets", {}).get("rows", [])
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
        # "disease" フィールドが null のエントリが先頭に来ることがあり（実際に
        # TP53 の TEPRASIRAN 等で確認）、フィルタせず diseases[0] を使うと
        # 実際は適応疾患があるのに空欄表示になってしまうため、null/空文字を除外する
        diseases = [
            (d.get("disease") or {}).get("name", "")
            for d in (row.get("diseases") or [])
            if isinstance(d, dict)
        ]
        diseases = [d for d in diseases if d]
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
        "disease_synonyms":   disease_synonyms,
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
