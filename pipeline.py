"""パイプライン — 1遺伝子の処理とバッチ実行をまとめる。

処理の流れ（遺伝子ごと）:
  1. collect_all        : 全DBから並列でエビデンス収集
  2. build_ppi_network  : IntAct + SIGNOR + Reactome から PPI ネットワーク構築
  3. run_network_enrichment : g:Profiler で機能エンリッチメント
  4. build_llm_context  : 収集結果を LLM 用コンテキストに整形
  5. generate_hypothesis: 仮説レポート本文を生成
  6. build_report       : Markdown レポートを組み立てて保存
"""
from __future__ import annotations

import os
import json
from datetime import datetime
from pathlib import Path

# BIOGRID_API_KEY が設定されていれば BioGRID を PPI に含める（非商用ライセンス）
_USE_BIOGRID = bool(os.environ.get("BIOGRID_API_KEY"))

import aggregator
import network as net
import hypothesis as hyp
import report
from collectors import uniprot

REPORTS_DIR = Path("reports")


def process_gene(
    gene: str,
    disease: str,
    disease_id: str,
    llm,
    lang: str = "en",
    context_config: dict | None = None,
    verbose: bool = True,
) -> dict:
    """1遺伝子×疾患を処理し、レポートを保存して結果サマリーを返す。

    Returns:
        {"gene", "status", "eval", "path"}（失敗時は path なし）
    """
    def log(msg=""):
        if verbose:
            print(msg)

    # 1. エビデンス収集
    try:
        evidence = aggregator.collect_all(gene, disease, verbose=verbose, disease_id=disease_id)
    except Exception as e:
        log(f"  ✗ データ収集失敗: {e}")
        return {"gene": gene, "status": f"データ収集失敗: {e}"}

    # 2-3. PPI ネットワーク + エンリッチメント
    ppi_graph, enrichment = _build_network(gene, log)

    # 4. LLM コンテキスト
    context = aggregator.build_llm_context(evidence, config=context_config)
    if ppi_graph:
        partners = net.rank_partners(ppi_graph, gene.upper())[:10]
        partner_fns = uniprot.get_functions_for_genes(partners) if partners else {}
        context += "\n\n" + net.network_summary_for_llm(
            ppi_graph, gene, enrichment, partner_functions=partner_fns)
    log(f"  コンテキスト: {len(context):,} 文字")

    # 5. 仮説生成（ストリーミング表示）
    log("  仮説生成中...\n")
    try:
        cb = (lambda tok: print(tok, end="", flush=True)) if verbose else None
        hypothesis = hyp.generate_hypothesis(gene, disease, context, llm, lang=lang, stream_callback=cb)
        log(f"\n  ✓ 仮説生成完了 ({len(hypothesis):,} 文字)")
    except Exception as e:
        log(f"\n  ✗ 仮説生成失敗: {e}")
        return {"gene": gene, "status": f"仮説生成失敗: {e}"}

    # LLM が自己流の "## References" を書いていれば除去し、エビデンス
    # コンテキストから確定的に生成したリンク付き References に差し替える。
    references_section = report.references_md(evidence.get("full_references") or {})
    if references_section:
        hypothesis = report.strip_llm_references(hypothesis) + "\n\n---\n\n" + references_section
        if verbose:
            print("\n\n---\n\n" + references_section)

    # 6. レポート保存
    path = _save_report(gene, disease, lang, hypothesis, context,
                        evidence, ppi_graph, enrichment, log)
    return {"gene": gene, "status": "✓ 完了", "path": str(path)}


def run_batch(
    genes: list[str],
    selected_disease: dict,
    llm,
    lang: str = "en",
    context_config: dict | None = None,
    verbose: bool = True,
) -> list[dict]:
    """複数遺伝子を順に処理し、結果サマリーのリストを返す。"""
    disease = selected_disease["name"]
    disease_id = selected_disease["id"]
    if verbose:
        print(f"疾患: {disease}  ({disease_id})")
        print(f"対象遺伝子 ({len(genes)}件): {', '.join(genes)}")
        print("=" * 60)

    results = []
    for i, gene in enumerate(genes, 1):
        if verbose:
            print(f"\n[{i}/{len(genes)}] {gene} × {disease}\n" + "-" * 50)
        results.append(process_gene(gene, disease, disease_id, llm, lang, context_config, verbose))
    return results


# ──────────────────────────────────────────────────────────────
# 内部ヘルパー
# ──────────────────────────────────────────────────────────────
def _build_network(gene: str, log):
    """PPI ネットワークとエンリッチメントを構築する（失敗しても None を返す）。"""
    log("  PPIネットワーク構築中...")
    try:
        graph = net.build_ppi_network(gene, use_biogrid=_USE_BIOGRID,
                                      use_reactome=False, use_intact=False)
        enrichment = net.run_network_enrichment(graph) if graph else {}
        if graph:
            log(f"  ✓ PPI: {graph.number_of_nodes()} nodes / {graph.number_of_edges()} edges")
            n_enr = len((enrichment or {}).get("results", []))
            if n_enr:
                log(f"  ✓ エンリッチメント: {n_enr} 有意項目")
        return graph, enrichment
    except Exception as e:
        log(f"  ⚠ ネットワーク構築エラー: {e}")
        return None, {}


def _save_report(gene, disease, lang, hypothesis, context,
                 evidence, ppi_graph, enrichment, log) -> Path:
    """レポート（.md）と付随データ（raw JSON, PPI画像）を保存し、md パスを返す。"""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    pair_dir = REPORTS_DIR / f"{gene}_{disease.replace(' ', '_')}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    # PPI 画像
    ppi_section = ""
    if ppi_graph and ppi_graph.number_of_edges() > 0:
        img_name = f"{ts}_ppi.png"
        partners = net.rank_partners(ppi_graph, gene.upper())[:30]
        if net.render_ppi_image(ppi_graph, gene, str(pair_dir / img_name),
                                enrichment=enrichment, max_nodes=30):
            ppi_section = report.ppi_md(gene, img_name, partners=partners)
            log(f"  ✓ PPI画像: {pair_dir / img_name}")

    # レポート本体
    md = report.build_report(
        gene, disease, lang, hypothesis, context,
        generated_iso=datetime.now().isoformat(),
        ppi_section=ppi_section,
        enrichment_section=report.enrichment_md(enrichment),
        competitive_section=report.competitive_landscape_md(
            evidence.get("evidence", {}).get("clinicaltrials")),
    )
    rpt_path = pair_dir / f"{ts}_{'JA' if lang == 'ja' else 'EN'}.md"
    rpt_path.write_text(md, encoding="utf-8")

    # 付随 JSON（生データ）
    (pair_dir / f"{ts}_raw.json").write_text(
        json.dumps(evidence, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    log(f"  ✓ 保存: {rpt_path}")
    return rpt_path


def save_summary(results: list[dict], disease: str) -> Path:
    """バッチサマリー（Markdown）を保存してパスを返す。"""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    md = report.summary_md(results, disease, datetime.now().strftime("%Y-%m-%d %H:%M"))
    path = REPORTS_DIR / f"{disease.replace(' ', '_')}_summary_{ts}.md"
    path.write_text(md, encoding="utf-8")
    return path
