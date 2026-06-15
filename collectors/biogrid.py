"""BioGRID PPI collector (MIT License / academic non-commercial use).

注意: BioGRID の利用規約は非商用・学術研究目的限定です。
APIキーは https://webservice.thebiogrid.org/ で無料登録取得できます。
"""
import os
import requests

BIOGRID_API = "https://webservice.thebiogrid.org/interactions/"

def get_interactions(gene_symbol: str, api_key: str = None) -> list[dict]:
    """Return BioGRID PPIs for gene_symbol.

    api_key が None の場合は環境変数 BIOGRID_API_KEY を参照。
    APIキーが設定されていなければ空リストを返し、警告を出す。
    """
    key = api_key or os.environ.get("BIOGRID_API_KEY", "")
    if not key:
        print("  [BioGRID] APIキー未設定 (BIOGRID_API_KEY)。スキップします。")
        print("  登録: https://webservice.thebiogrid.org/")
        return []

    params = {
        "accessKey":          key,
        "geneList":           gene_symbol,
        "searchNames":        "true",
        "includeHeader":      "true",
        "taxId":              "9606",
        "interSpeciesExcluded": "true",
        "selfInteractionsExcluded": "true",
        "format":             "json",
        "max":                200,
        "start":              0,
    }

    try:
        r = requests.get(BIOGRID_API, params=params, timeout=20)
        r.raise_for_status()
        data = r.json()
    except requests.HTTPError as e:
        if r.status_code == 400:
            print(f"  [BioGRID] 400エラー: APIキーが正しくない可能性があります: {r.text[:200]}")
        else:
            print(f"  [BioGRID] HTTPエラー: {e}")
        return []
    except Exception as e:
        print(f"  [BioGRID] エラー: {e}")
        return []

    results = []
    gene_upper = gene_symbol.upper()

    for interaction_id, item in data.items():
        sym_a = item.get("OFFICIAL_SYMBOL_A", "")
        sym_b = item.get("OFFICIAL_SYMBOL_B", "")
        partner = sym_b if sym_a.upper() == gene_upper else sym_a
        exp_system = item.get("EXPERIMENTAL_SYSTEM", "")
        exp_type   = item.get("EXPERIMENTAL_SYSTEM_TYPE", "")
        pubmed_id  = str(item.get("PUBMED_ID", ""))
        score      = item.get("SCORE", None)

        results.append({
            "source":      sym_a,
            "target":      sym_b,
            "partner":     partner,
            "direction":   "—",  # BioGRID は無向グラフ
            "effect":      "physical association",
            "mechanism":   exp_system,
            "exp_type":    exp_type,
            "pmid":        pubmed_id,
            "score":       float(score) if score not in (None, "") else None,
            "db":          "BioGRID",
        })

    return results
