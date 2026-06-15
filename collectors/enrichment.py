"""GO/Pathway エンリッチメント解析 via g:Profiler REST API (BSD 2-Clause).

APIキー不要。商用利用可。
エンドポイント: https://biit.cs.ut.ee/gprofiler/api/gost/profile/
ソース: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP
"""
import requests

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

SOURCES = ["GO:BP", "GO:MF", "GO:CC", "REAC", "WP"]


def run_enrichment(
    gene_list: list[str],
    organism: str = "hsapiens",
    top_n: int = 30,
) -> list[dict]:
    """g:Profiler でエンリッチメント解析を実行し、上位 top_n 件を返す。

    Args:
        gene_list: HGNC シンボルのリスト
        organism:  g:Profiler 生物種 ID (デフォルト: hsapiens)
        top_n:     返す結果の最大件数

    Returns:
        [{source, term_id, term_name, p_value, intersection_size, gene_ratio, genes}]
    """
    if not gene_list:
        return []

    payload = {
        "organism":        organism,
        "query":           gene_list,
        "sources":         SOURCES,
        "user_threshold":  0.05,
        "significance_threshold_method": "fdr",
        "no_evidences":    True,
    }

    try:
        r = requests.post(GPROFILER_URL, json=payload, timeout=30)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"  [Enrichment] g:Profiler エラー: {e}")
        return []

    results_raw = (data.get("result") or [])
    if not results_raw:
        return []

    results = []
    for item in results_raw[:top_n]:
        query_size       = item.get("query_size", 1) or 1
        intersection_size = item.get("intersection_size", 0)
        results.append({
            "source":           item.get("source", ""),
            "term_id":          item.get("native", ""),
            "term_name":        item.get("name", ""),
            "p_value":          item.get("p_value", 1.0),
            "intersection_size": intersection_size,
            "gene_ratio":       intersection_size / query_size,
            "genes":            item.get("intersections", []),
        })

    # p_value 昇順でソート
    results.sort(key=lambda x: x["p_value"])

    return results


def top_terms_by_source(enrichment_results: list[dict], top_per_source: int = 5) -> dict:
    """ソース別に上位 top_per_source 件を返す辞書。"""
    grouped: dict[str, list] = {}
    for item in enrichment_results:
        src = item["source"]
        if src not in grouped:
            grouped[src] = []
        if len(grouped[src]) < top_per_source:
            grouped[src].append(item)
    return grouped
