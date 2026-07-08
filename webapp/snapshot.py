"""Web アプリの表示内容をオフライン閲覧可能な単一 HTML として書き出す。

index.html の描画ロジック（JS の createGeneCard / updateCollectorChip /
handleEvent など）をそのまま再利用するため、テンプレートに
<script id="snapshot-data" type="application/json"> でデータを埋め込むだけの
軽量な仕組みにしている。生成された HTML は、サーバーを介さずブラウザで直接
開いて「取得データ／仮説／PPI／エンリッチメント／ダウンロード」の全タブを
その場で閲覧できる（レポート済みの静的スナップショット）。
"""
import json
from pathlib import Path

_TEMPLATE_PATH = Path(__file__).parent / "templates" / "index.html"


def build_snapshot_html(
    gene: str,
    disease: str,
    lang: str,
    generated_iso: str,
    collectors: dict,
    hypothesis: str,
    ppi_image_filename: str,
    partners: list[str],
    partner_functions: list[dict],
    enrichment_results: list[dict],
    excluded_hubs: list[dict],
    report_filename: str,
) -> str:
    """スナップショット HTML 文字列を生成する。

    collectors: {source: {"ok": bool, "summary": str, "data": dict|None}}
        取得データタブの各コレクター結果（webapp.app の results/errors から整形）。
    ppi_image_filename / report_filename: スナップショット HTML と同じ
        ディレクトリに保存されているファイル名（相対パスで参照するため）。
    """
    template = _TEMPLATE_PATH.read_text(encoding="utf-8")

    snapshot = {
        "gene": gene,
        "disease": disease,
        "lang": lang,
        "generated": generated_iso,
        "collectors": collectors,
        "hypothesis": hypothesis,
        "ppi_image": ppi_image_filename,
        "partners": partners,
        "partner_functions": partner_functions,
        "enrichment_results": enrichment_results,
        "excluded_hubs": excluded_hubs,
        "report_filename": report_filename,
    }
    payload = json.dumps(snapshot, ensure_ascii=False, default=str)
    script_tag = f'<script id="snapshot-data" type="application/json">{payload}</script>\n'

    if "</body>" in template:
        html = template.replace("</body>", script_tag + "</body>")
    else:
        html = template + script_tag
    return html
