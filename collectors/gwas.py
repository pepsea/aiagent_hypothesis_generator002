"""GWAS Catalog REST API (EMBL-EBI, free for all use)."""
import requests
from concurrent.futures import ThreadPoolExecutor

BASE = "https://www.ebi.ac.uk/gwas/rest/api"


def _fetch_snp_hits(snp: dict) -> list[dict]:
    """1 SNP の associations と trait 名を取得（並列実行用）。"""
    rsid = snp.get("rsId", "")
    assoc_link = (snp.get("_links", {}).get("associations") or {}).get("href")
    if not assoc_link:
        return []
    try:
        ra = requests.get(assoc_link, timeout=12)
        ra.raise_for_status()
        associations = ra.json().get("_embedded", {}).get("associations", [])
    except Exception:
        return []

    hits = []
    for assoc in associations:
        trait = ""
        gwas_url = ""
        tl = (assoc.get("_links", {}).get("efoTraits") or {}).get("href")
        if tl:
            try:
                rt = requests.get(tl, timeout=8)
                traits = rt.json().get("_embedded", {}).get("efoTraits", [])
                trait = ", ".join(t.get("trait", "") for t in traits if t.get("trait"))
                # GWAS Catalog の trait ページ（当該形質の全アソシエーション一覧）にリンク
                first_efo = next((t.get("shortForm") for t in traits if t.get("shortForm")), None)
                if first_efo:
                    gwas_url = f"https://www.ebi.ac.uk/gwas/efotraits/{first_efo}"
            except Exception:
                pass
        hits.append({
            "trait":                 trait,
            "p_value":               assoc.get("pvalue"),
            "or_beta":               assoc.get("orPerCopyNum") or assoc.get("betaNum"),
            "snps":                  [rsid],
            "risk_allele_frequency": assoc.get("riskFrequency"),
            "gwas_url":              gwas_url,
        })
    return hits


def get_gwas_associations(gene_symbol: str, disease_query: str = None,
                          max_snps: int = 15, max_results: int = 20) -> list[dict]:
    """Return GWAS hits for a gene, optionally filtered by trait.

    旧 /genes/{gene}/associations は廃止 (500) されたため、
    findByGene で SNP を取得し、各 SNP の associations → efoTraits を並列で辿る。
    """
    try:
        r = requests.get(
            f"{BASE}/singleNucleotidePolymorphisms/search/findByGene",
            params={"geneName": gene_symbol, "size": max_snps}, timeout=20)
        if r.status_code == 404:
            return []
        r.raise_for_status()
        snps = r.json().get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
    except Exception:
        return []

    # SNP ごとの取得を並列化（逐次だと ~60s → 並列で数秒）
    # 語順違い（例: "type 2 diabetes mellitus" vs "diabetes mellitus type 2"）で
    # 単純な部分文字列一致だと本来ヒットすべき trait を取りこぼすため、
    # クエリを単語分割し全単語がトレイト文字列に含まれるかで判定する
    query_words = disease_query.lower().split() if disease_query else []

    results, seen = [], set()
    with ThreadPoolExecutor(max_workers=10) as ex:
        for hits in ex.map(_fetch_snp_hits, snps):
            for h in hits:
                trait_lower = (h["trait"] or "").lower()
                if query_words and not all(w in trait_lower for w in query_words):
                    continue
                key = (h["snps"][0], h["trait"])
                if key in seen:
                    continue
                seen.add(key)
                results.append(h)

    def _pv(x):
        try:
            return float(x.get("p_value") or 1)
        except (TypeError, ValueError):
            return 1.0
    results.sort(key=_pv)
    return results[:max_results]


# ClinVar の trait_name にしばしば入るプレースホルダー（実際の疾患名ではない）
_CLINVAR_NO_CONDITION = {"", "not provided", "not specified", "see cases"}


def get_clinvar_variants(
    gene_symbol: str,
    disease_query: str = None,
    disease_synonyms: list = None,
    max_results: int = 100,
) -> list[dict]:
    """Return ClinVar pathogenic variants for a gene, filtered by disease.

    疾患名が分類されている（trait_name が実際の疾患名であり、"not provided"等の
    プレースホルダーではない）レコードを中心に抽出する。候補を多めに取得して
    フィルタしたうえで、最終評価日（last_evaluated）降順で直近優先に並べ、
    最大 max_results 件（デフォルト100）を返す。

    disease_query: 入力疾患名。指定時は esearch クエリと condition ポストフィルタに使用。
    disease_synonyms: OpenTargets から取得した疾患の同義語リスト（表記ゆらぎ対応）。
    """
    import time
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # 疾患名を esearch クエリに含めて API 側で絞り込む
    disease_term = f' AND "{disease_query}"[Disease/Phenotype]' if disease_query else ""
    query = (
        f"{gene_symbol}[Gene Name]"
        f" AND (Pathogenic[Clinical significance] OR Likely pathogenic[Clinical significance])"
        f"{disease_term}"
    )

    # フィルタで減る分の余裕を持たせて多めに取得
    r = requests.get(f"{base}/esearch.fcgi", params={
        "db": "clinvar", "term": query, "retmax": max_results * 3, "retmode": "json"
    }, timeout=15)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])

    if not ids:
        return []

    # 429対策: リトライ付きで esummary を呼ぶ（POST: ID数が多い場合のURL長対策）
    for attempt in range(3):
        time.sleep(1 + attempt * 2)  # 1s, 3s, 5s
        r2 = requests.post(f"{base}/esummary.fcgi", data={
            "db": "clinvar", "id": ",".join(ids), "retmode": "json"
        }, timeout=20)
        if r2.status_code == 429:
            continue
        r2.raise_for_status()
        break
    else:
        return []  # リトライ上限に達した場合は空を返す
    result = r2.json().get("result", {})

    # ポストフィルタ用キーワードセット（入力名 + synonyms + 入力名の個別トークン）
    filter_terms = []
    if disease_query:
        filter_terms.append(disease_query.lower())
        filter_terms.extend(t.lower() for t in disease_query.split() if len(t) > 3)
    for syn in (disease_synonyms or []):
        filter_terms.append(syn.lower())

    variants = []
    for vid in ids:
        if vid not in result:
            continue
        item = result[vid]
        # NCBI eutils は germline_classification（旧 clinical_significance）に
        # 臨床的意義・レビュー状態・trait を格納する
        cls = item.get("germline_classification") or item.get("clinical_significance") or {}
        traits = cls.get("trait_set") or []
        condition = traits[0].get("trait_name", "") if traits else ""
        if condition.strip().lower() in _CLINVAR_NO_CONDITION:
            continue  # 疾患名が分類されていないレコードは対象外
        # 疾患フィルタ: condition または title に疾患名・synonym が含まれるか確認
        if filter_terms:
            target_text = (condition + " " + item.get("title", "")).lower()
            if not any(term in target_text for term in filter_terms):
                continue
        variants.append({
            "variant_id": vid,
            "title": item.get("title", ""),
            "clinical_significance": cls.get("description", ""),
            "condition": condition,
            "review_status": cls.get("review_status", ""),
            "last_evaluated": cls.get("last_evaluated", ""),
        })

    variants.sort(key=lambda v: v.get("last_evaluated") or "", reverse=True)
    return variants[:max_results]
