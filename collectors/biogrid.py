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

    # BioGRID は同一ペアについて実験手法や文献違いで複数レコードを返すことが多い
    # （異なる EXPERIMENTAL_SYSTEM / PUBMED_ID）。取得データ表示・PPIネットワーク
    # 構築の両方で遺伝子ごとに重複させないため、パートナーごとに集約し1件のみ残す。
    # 選択基準: スコアがあれば最大のもの、無ければ最初に見つかった1件。
    best_by_partner: dict[str, dict] = {}
    gene_upper = gene_symbol.upper()

    for interaction_id, item in data.items():
        sym_a = item.get("OFFICIAL_SYMBOL_A", "")
        sym_b = item.get("OFFICIAL_SYMBOL_B", "")
        partner = sym_b if sym_a.upper() == gene_upper else sym_a
        if not partner or partner.upper() == gene_upper:
            continue

        exp_system = item.get("EXPERIMENTAL_SYSTEM", "")
        exp_type   = item.get("EXPERIMENTAL_SYSTEM_TYPE", "")
        pubmed_id  = str(item.get("PUBMED_ID", ""))
        score_raw  = item.get("SCORE", None)
        score      = float(score_raw) if score_raw not in (None, "") else None

        key = partner.upper()
        existing = best_by_partner.get(key)
        if existing is not None:
            existing_score = existing.get("score")
            # 既存にスコアがあり、新規のスコアがそれ以下なら置き換えない
            if existing_score is not None and (score is None or score <= existing_score):
                continue
            # 既存にスコアが無く、新規にも無い場合は既存（先着）を維持
            if existing_score is None and score is None:
                continue

        best_by_partner[key] = {
            "source":      sym_a,
            "target":      sym_b,
            "partner":     partner,
            "direction":   "—",  # BioGRID は無向グラフ
            "effect":      "physical association",
            "mechanism":   exp_system,
            "exp_type":    exp_type,
            "pmid":        pubmed_id,
            "score":       score,
            "db":          "BioGRID",
        }

    results = sorted(best_by_partner.values(),
                     key=lambda x: x.get("score") if x.get("score") is not None else -1,
                     reverse=True)
    return results
