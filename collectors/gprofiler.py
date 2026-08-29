"""g:Profiler — pathway enrichment for a gene list (free, no key required).

POST https://biit.cs.ut.ee/gprofiler/api/gost/profile/
Returns enriched terms from Reactome (REAC), KEGG, and GO:BP.
"""
import requests

GPROFILER_API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

def enrich_gene_list(
    gene_symbols: list[str],
    organism: str = "hsapiens",
    sources: list[str] = None,
    significance_threshold: float = 0.05,
    max_results: int = 30,
) -> list[dict]:
    """Run g:Profiler enrichment on a list of gene symbols.

    Returns list of enriched pathways sorted by p-value:
    [{"source": "REAC", "term_id": "R-HSA-...", "name": "...", "p_value": 0.001,
      "intersection_size": 5, "term_size": 120, "genes": [...]}]
    """
    if not gene_symbols:
        return []
    # KEGG は商用利用不可のため除外。REAC / GO:BP / WP のみ使用。
    if sources is None:
        sources = ["REAC", "GO:BP", "WP"]

    try:
        payload = {
            "organism": organism,
            "query": gene_symbols,
            "sources": sources,
            "user_threshold": significance_threshold,
            "ordered": False,
            "all_results": False,
        }
        r = requests.post(GPROFILER_API, json=payload, timeout=30)
        r.raise_for_status()
        data = r.json()
    except Exception:
        return []

    results_raw = data.get("result") or []
    meta = data.get("meta") or {}
    # クエリ遺伝子の順序リスト（intersections のインデックス対応）
    queries = (meta.get("query_metadata") or {}).get("queries") or {}
    ordered_genes = queries.get("query_1", gene_symbols)

    results = []
    for item in results_raw:
        # intersections[i] が非空ならi番目のクエリ遺伝子がヒット
        raw_ix = item.get("intersections") or []
        hit_genes = [
            ordered_genes[i]
            for i, ann in enumerate(raw_ix)
            if ann and i < len(ordered_genes)
        ]
        results.append({
            "source":            item.get("source", ""),
            "term_id":           item.get("native", ""),
            "name":              item.get("name", ""),
            "p_value":           item.get("p_value", 1.0),
            "intersection_size": item.get("intersection_size", 0),
            "term_size":         item.get("term_size", 0),
            "genes":             hit_genes,
        })

    results.sort(key=lambda x: x["p_value"])
    return results[:max_results]
