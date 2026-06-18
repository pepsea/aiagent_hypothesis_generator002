"""レポート生成 — 仮説・評価・PPI・エンリッチメントを Markdown / HTML に整形する。

このモジュールは「文字列を組み立てて返す」だけに徹し、ファイル保存や画面表示
（display）は呼び出し側（pipeline / notebook）が行う。
"""
from __future__ import annotations

import re
from collections import defaultdict

# ──────────────────────────────────────────────────────────────
# Target Validity 評価（4段階）
# ──────────────────────────────────────────────────────────────
# 評価項目（キー, 表示ラベル）— hypothesis.generate_presentation_eval の出力に対応
EVAL_KEYS = [
    ("genetic_association",     "Genetic Association"),
    ("functional_association",  "Functional Association"),
    ("clinical_relevance",      "Clinical Relevance"),
    ("network_context",         "Expression / Network"),
    ("target_validity_overall", "Overall"),
]

# 4段階レーティングの略号と色
MARK_COLOR = {"VH": "#0550ae", "H": "#1a7f37", "M": "#9a6700", "L": "#cf222e", "—": "#aaa"}
MARK_BG    = {"VH": "#ddf4ff", "H": "#e6ffed", "M": "#fff8c5", "L": "#ffebe9", "—": "#f6f8fa"}

_PLACEHOLDER_FINDINGS = {
    "...", "FINDING", "finding", "RATING",
    "一文", "総合評価を一文で", "仮説の一文要約",
    "one sentence hypothesis summary", "one-sentence overall summary",
}


def to_mark(rating: str) -> str:
    """rating 文字列を VH / H / M / L / — に変換する。"""
    s = (rating or "").strip().lower()
    if not s or s in ("...", "-", "rating"):
        return "—"
    if re.search(r"very.?high|非常に高", s):     return "VH"
    if re.search(r"\bhigh\b|^高", s):            return "H"
    if re.search(r"middle|moderate|中", s):      return "M"
    if re.search(r"\blow\b|^低|weak|弱", s):     return "L"
    return "—"


def get_mark(ev: dict, key: str) -> tuple[str, str]:
    """評価 dict から (略号, finding) を取り出す。"""
    item    = (ev or {}).get(key) or {}
    mark    = to_mark(item.get("rating", ""))
    finding = (item.get("finding") or "").strip()
    if finding in _PLACEHOLDER_FINDINGS:
        finding = ""
    return mark, finding


_VALIDITY_LEGEND = (
    "<sub><sup>VH=Very High (clear association &amp; severity-linked) &nbsp;|&nbsp; "
    "H=High (clear association, pathway-direct) &nbsp;|&nbsp; "
    "M=Middle (explainable via PPI/pathway) &nbsp;|&nbsp; "
    "L=Low (no data)</sup></sub>"
)


def validity_table_md(ev: dict) -> str:
    """個別レポート冒頭の Target Validity テーブル（Markdown）。"""
    lines = ["## Target Validity", "", "| Evaluation | Rating | Finding |", "|---|:---:|---|"]
    for key, label in EVAL_KEYS:
        mark, finding = get_mark(ev, key)
        b = "**" if key == "target_validity_overall" else ""
        lines.append(f"| {b}{label}{b} | {b}{mark}{b} | {finding} |")
    lines += ["", _VALIDITY_LEGEND, ""]
    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────
# 機能エンリッチメント
# ──────────────────────────────────────────────────────────────
# g:Profiler ソース ID → 表示名（KEGG/TF など商用ソースは使用しない）
_SOURCE_LABEL = {
    "GO:BP":  "GO Biological Process",
    "GO:MF":  "GO Molecular Function",
    "GO:CC":  "GO Cellular Component",
    "REAC":   "Reactome",
    "WP":     "WikiPathways",
    "HP":     "Human Phenotype",
    "CORUM":  "Protein Complexes",
    "HPA":    "Human Protein Atlas",
}
_SOURCE_ORDER = ["GO:BP", "GO:MF", "GO:CC", "REAC", "WP", "HP", "CORUM", "HPA"]


def enrichment_md(enrichment: dict, top_per_source: int = 5) -> str:
    """エンリッチメント結果をソース別テーブルにまとめた Markdown。結果なしなら空文字。"""
    results = (enrichment or {}).get("results", [])
    if not results:
        return ""

    by_source: dict[str, list] = defaultdict(list)
    for r in results:
        by_source[r["source"]].append(r)

    ordered = [(s, by_source[s]) for s in _SOURCE_ORDER if s in by_source]
    ordered += [(s, by_source[s]) for s in sorted(by_source) if s not in _SOURCE_ORDER]

    lines = ["## Functional Enrichment (g:Profiler, FDR < 0.05)", ""]
    for source, terms in ordered:
        lines.append(f"### {_SOURCE_LABEL.get(source, source)}")
        lines.append("")
        lines.append("| Term | p-value | Genes (overlap) |")
        lines.append("|---|:---:|---|")
        for t in terms[:top_per_source]:
            genes = ", ".join(t.get("genes", [])[:8])
            lines.append(
                f"| {t['term_name'][:70]} `{t.get('term_id', '')}` "
                f"| {t['p_value']:.2e} | {genes} ({t.get('intersection_size', 0)}) |"
            )
        lines.append("")
    lines.append(
        f"<sub><sup>有意項目合計: {len(results)} 件 (FDR&lt;0.05) — "
        f"[g:Profiler](https://biit.cs.ut.ee/gprofiler/gost)</sup></sub>\n"
    )
    return "\n".join(lines)


def ppi_md(gene: str, image_filename: str) -> str:
    """PPI ネットワーク画像を埋め込む Markdown セクション。"""
    return (
        f"## PPI Network\n\n"
        f"![PPI network of {gene}]({image_filename})\n\n"
        f"<sub><sup>★ = {gene} (target, 上部) ／ 下部 = PPI パートナー"
        f"（色 = enrichment 上位パスウェイ）</sup></sub>\n\n"
    )


# ──────────────────────────────────────────────────────────────
# 個別レポート全体
# ──────────────────────────────────────────────────────────────
def build_report(
    gene: str,
    disease: str,
    lang: str,
    eval_result: dict,
    hypothesis: str,
    context: str,
    generated_iso: str,
    ppi_section: str = "",
    enrichment_section: str = "",
) -> str:
    """1遺伝子×疾患の完全な Markdown レポートを組み立てて返す。"""
    parts = [
        f"# Drug Discovery Hypothesis: {gene} × {disease}",
        f"Generated: {generated_iso}  |  Language: {lang}",
        "", "---", "",
        validity_table_md(eval_result),
        "---", "",
    ]
    if ppi_section:
        parts.append(ppi_section)
    if enrichment_section:
        parts += [enrichment_section, "---", ""]
    parts += [hypothesis, "", "---", "", "## Evidence Context", "", context]
    return "\n".join(parts)


# ──────────────────────────────────────────────────────────────
# バッチサマリー（行=評価項目、列=遺伝子）
# ──────────────────────────────────────────────────────────────
def summary_html(results: list[dict], disease: str) -> str:
    """評価サマリーを HTML テーブルで返す（Jupyter 表示用）。"""
    legend = (
        '<div style="display:flex;gap:14px;flex-wrap:wrap;font-size:12px;color:#444;margin:0 0 10px">'
        '<span><b style="background:#ddf4ff;color:#0550ae;padding:2px 7px;border-radius:3px">VH</b> Very High</span>'
        '<span><b style="background:#e6ffed;color:#1a7f37;padding:2px 7px;border-radius:3px">H</b> High</span>'
        '<span><b style="background:#fff8c5;color:#9a6700;padding:2px 7px;border-radius:3px">M</b> Middle</span>'
        '<span><b style="background:#ffebe9;color:#cf222e;padding:2px 7px;border-radius:3px">L</b> Low</span>'
        '</div>'
    )
    gene_th = "".join(
        f'<th style="padding:6px 16px;text-align:center;border:1px solid #ddd;background:#f5f5f5">{r["gene"]}</th>'
        for r in results
    )
    rows = ""
    for key, label in EVAL_KEYS:
        overall = key == "target_validity_overall"
        row_bg  = "background:#ececec;" if overall else ""
        lbl_st  = "font-weight:bold;" if overall else ""
        cells = ""
        for r in results:
            mark, tip = get_mark(r.get("eval", {}), key)
            title = f' title="{tip}"' if tip else ""
            cells += (
                f'<td style="padding:7px 16px;text-align:center;border:1px solid #ddd;'
                f'background:{MARK_BG.get(mark, "#f6f8fa")};font-weight:bold;'
                f'color:{MARK_COLOR.get(mark, "#aaa")};font-size:15px"{title}>{mark}</td>'
            )
        rows += (
            f'<tr style="{row_bg}"><td style="padding:7px 12px;border:1px solid #ddd;'
            f'white-space:nowrap;{lbl_st}">{label}</td>{cells}</tr>'
        )
    return (
        f'<h3 style="margin-bottom:6px">Target Validity Summary — {disease}</h3>{legend}'
        f'<table style="border-collapse:collapse;font-size:13px"><thead><tr>'
        f'<th style="padding:6px 12px;text-align:left;border:1px solid #ddd;background:#f5f5f5">Evaluation</th>'
        f'{gene_th}</tr></thead><tbody>{rows}</tbody></table>'
        f'<p style="font-size:11px;color:#888;margin:6px 0 0">※ セルにカーソルを合わせると根拠が表示されます</p>'
    )


def summary_md(results: list[dict], disease: str, generated: str) -> str:
    """評価サマリーを Markdown で返す（ファイル保存用）。"""
    genes = [r["gene"] for r in results]
    lines = [
        f"# Target Validity Summary — {disease}",
        f"Generated: {generated}", "",
        "VH=Very High  H=High  M=Middle  L=Low  —=No data", "",
        "| Evaluation | " + " | ".join(genes) + " |",
        "|------------|" + "|".join([":----------:"] * len(results)) + "|",
    ]
    for key, label in EVAL_KEYS:
        marks = [get_mark(r.get("eval", {}), key)[0] for r in results]
        lines.append(f"| **{label}** | " + " | ".join(marks) + " |")
    return "\n".join(lines)
