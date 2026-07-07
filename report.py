"""レポート生成 — 仮説・評価・PPI・エンリッチメントを Markdown / HTML に整形する。

このモジュールは「文字列を組み立てて返す」だけに徹し、ファイル保存や画面表示
（display）は呼び出し側（pipeline / notebook）が行う。
"""
from __future__ import annotations

from collections import defaultdict

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


def ppi_md(gene: str, image_filename: str, partners: list[str] | None = None,
           functions: dict | None = None) -> str:
    """PPI ネットワーク画像 + パートナー遺伝子リスト + 機能情報の Markdown セクション。

    functions: {GENE(upper): {"protein_name","function"}} を渡すと、
        対象遺伝子および PPI パートナーの UniProt 機能説明をリスト表示する。
    """
    lines = [
        "## PPI Network",
        "",
        f"![PPI network of {gene}]({image_filename})",
        "",
        f"<sub><sup>★ = {gene} (target, 上部) ／ 下部 = PPI パートナー。"
        f"色 = データソース（SIGNOR / BioGRID、複数DB共通は強調色）"
        f"</sup></sub>",
        "",
    ]
    if partners:
        lines += [
            f"**PPI Partners ({len(partners)} genes):** "
            + ", ".join(f"`{p}`" for p in partners),
            "",
        ]

    # 対象遺伝子・PPI 遺伝子の機能情報（UniProt）
    if functions:
        lines += ["### タンパク質機能 (UniProt)", ""]
        order = [gene.upper()] + [p for p in (partners or []) if p.upper() != gene.upper()]
        seen = set()
        for g in order:
            info = functions.get(g.upper())
            if not info or g.upper() in seen:
                continue
            seen.add(g.upper())
            func = (info.get("function") or "").strip()
            if not func:
                continue
            pname = info.get("protein_name", "")
            label = f"{g}" + (f" — {pname}" if pname else "")
            tag = " **(target)**" if g.upper() == gene.upper() else ""
            lines.append(f"- **{label}**{tag}: {func}")
        lines.append("")

    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────
# 個別レポート全体
# ──────────────────────────────────────────────────────────────
def build_report(
    gene: str,
    disease: str,
    lang: str,
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
        hypothesis, "",
    ]
    if ppi_section or enrichment_section:
        parts += ["---", "", "## Supporting Evidence", ""]
        if ppi_section:
            parts.append(ppi_section)
        if enrichment_section:
            parts.append(enrichment_section)
    parts += ["---", "", "## Evidence Context", "", context]
    return "\n".join(parts)


# ──────────────────────────────────────────────────────────────
# バッチサマリー（遺伝子ごとの処理結果一覧）
# ──────────────────────────────────────────────────────────────
def summary_html(results: list[dict], disease: str) -> str:
    """バッチ処理結果を HTML テーブルで返す（Jupyter 表示用）。"""
    rows = ""
    for r in results:
        done = r.get("status", "").startswith("✓")
        color = "#1a7f37" if done else "#cf222e"
        path  = r.get("path", "")
        rows += (
            f'<tr>'
            f'<td style="padding:6px 14px;border:1px solid #ddd;font-weight:bold">{r["gene"]}</td>'
            f'<td style="padding:6px 14px;border:1px solid #ddd;color:{color}">{r.get("status","")}</td>'
            f'<td style="padding:6px 14px;border:1px solid #ddd;color:#555;font-size:12px">{path}</td>'
            f'</tr>'
        )
    n_done = sum(1 for r in results if r.get("status", "").startswith("✓"))
    return (
        f'<h3 style="margin-bottom:6px">Batch Summary — {disease}</h3>'
        f'<p style="font-size:13px;color:#444;margin:0 0 8px">完了: {n_done}/{len(results)} 遺伝子</p>'
        f'<table style="border-collapse:collapse;font-size:13px"><thead><tr>'
        f'<th style="padding:6px 14px;text-align:left;border:1px solid #ddd;background:#f5f5f5">遺伝子</th>'
        f'<th style="padding:6px 14px;text-align:left;border:1px solid #ddd;background:#f5f5f5">ステータス</th>'
        f'<th style="padding:6px 14px;text-align:left;border:1px solid #ddd;background:#f5f5f5">レポート</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>'
    )


def summary_md(results: list[dict], disease: str, generated: str) -> str:
    """バッチ処理結果を Markdown で返す（ファイル保存用）。"""
    n_done = sum(1 for r in results if r.get("status", "").startswith("✓"))
    lines = [
        f"# Batch Summary — {disease}",
        f"Generated: {generated}", "",
        f"完了: {n_done}/{len(results)} 遺伝子", "",
        "| 遺伝子 | ステータス | レポート |",
        "|---|---|---|",
    ]
    for r in results:
        lines.append(f"| {r['gene']} | {r.get('status','')} | {r.get('path','')} |")
    return "\n".join(lines)
