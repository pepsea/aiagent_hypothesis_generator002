"""GO/Pathway エンリッチメント解析 via g:Profiler REST API (BSD 2-Clause).

APIキー不要。商用利用可。
エンドポイント: https://biit.cs.ut.ee/gprofiler/api/gost/profile/
ソース: GO:BP, GO:MF, GO:CC, REAC, WP
  ※ KEGG・TRANSFAC(TF) 等の商用ライセンスが必要なソースは使用しない。

API 仕様変更 (2025):
  - result は flat list（ネスト無し）
  - intersections[i] は i番目クエリ遺伝子のソースアノテーションリスト
  - ヒット遺伝子名は meta.query_metadata.queries.query_1 との対応で取得
  - no_evidences: True を使うと intersections が返らないため削除
"""
import requests

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

SOURCES = ["GO:BP", "GO:MF", "GO:CC", "REAC", "WP"]


def run_enrichment(
    gene_list: list[str],
    organism: str = "hsapiens",
    top_n: int = 50,
) -> list[dict]:
    """g:Profiler でエンリッチメント解析を実行し、上位 top_n 件を返す。

    Returns:
        [{source, term_id, term_name, p_value, term_size,
          intersection_size, gene_ratio, genes}]
          genes: ヒットした遺伝子シンボルのリスト
    """
    if not gene_list:
        return []

    payload = {
        "organism":  organism,
        "query":     gene_list,
        "sources":   SOURCES,
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        # no_evidences を除去することで intersections（ヒット遺伝子情報）を取得
    }

    try:
        r = requests.post(GPROFILER_URL, json=payload, timeout=30)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"  [Enrichment] g:Profiler エラー: {e}")
        return []

    results_raw = data.get("result") or []
    meta        = data.get("meta") or {}

    # クエリ遺伝子の順序リスト（intersections のインデックス対応）
    queries = (meta.get("query_metadata") or {}).get("queries") or {}
    ordered_genes = queries.get("query_1", gene_list)

    # 有意な結果のみ
    sig = [item for item in results_raw if item.get("significant")]
    sig.sort(key=lambda x: x.get("p_value", 1.0))

    results = []
    for item in sig[:top_n]:
        query_size        = item.get("query_size", 1) or 1
        intersection_size = item.get("intersection_size", 0)

        # intersections[i] は非空リスト → i番目の遺伝子がヒット
        raw_ix = item.get("intersections") or []
        hit_genes = [
            ordered_genes[i]
            for i, ann in enumerate(raw_ix)
            if ann and i < len(ordered_genes)
        ]

        results.append({
            "source":           item.get("source", ""),
            "term_id":          item.get("native", ""),
            "term_name":        item.get("name", ""),
            "p_value":          item.get("p_value", 1.0),
            "term_size":        item.get("term_size", 0),
            "intersection_size": intersection_size,
            "gene_ratio":       intersection_size / query_size,
            "genes":            hit_genes,
        })

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
