"""PPI ネットワーク構築 + エンリッチメント解析統合モジュール.

データソース: IntAct (EMBL-EBI) + SIGNOR (UNIROMA2) + BioGRID (optional)
エンリッチメント: g:Profiler (GO/KEGG/Reactome/WikiPathways)
可視化: pyvis (インタラクティブHTML) または networkx/matplotlib
"""
from __future__ import annotations

import os
from typing import Optional

try:
    import networkx as nx
    HAS_NX = True
except ImportError:
    HAS_NX = False

from collectors import intact, signor, biogrid, enrichment as enrich_mod


# ──────────────────────────────────────────────────────────────
# 1. ネットワーク構築
# ──────────────────────────────────────────────────────────────

def build_ppi_network(
    gene_symbol: str,
    use_biogrid: bool = True,
    biogrid_api_key: str | None = None,
) -> Optional["nx.Graph"]:
    """IntAct + SIGNOR (+ BioGRID) の相互作用から NetworkX グラフを構築する。

    Returns:
        nx.Graph: ノード属性に db, color, direct_partner を持つグラフ
                  networkx 未インストールの場合は None
    """
    if not HAS_NX:
        print("  [Network] networkx が未インストールです: pip install networkx")
        return None

    G = nx.Graph()
    center = gene_symbol.upper()
    G.add_node(center, color="#FF6B6B", size=25, db="center", direct_partner=True)

    def add_edges(interactions: list[dict], source_label: str):
        for item in interactions:
            src = (item.get("source") or "").strip().upper()
            tgt = (item.get("target") or "").strip().upper()
            if not src or not tgt or src == tgt:
                continue
            partner = tgt if src == center else src
            if not partner:
                continue

            # ノード追加
            if partner not in G:
                G.add_node(partner, color=_db_color(source_label),
                           size=15, db=source_label, direct_partner=True)
            elif source_label not in G.nodes[partner].get("db", ""):
                G.nodes[partner]["db"] += f",{source_label}"

            # エッジ追加 or 属性更新
            if G.has_edge(center, partner):
                ed = G.edges[center, partner]
                ed["weight"] = ed.get("weight", 1) + 1
                if source_label not in ed.get("db", ""):
                    ed["db"] += f",{source_label}"
            else:
                G.add_edge(center, partner,
                           weight=1,
                           effect=item.get("effect", ""),
                           mechanism=item.get("mechanism", ""),
                           db=source_label)

    # --- IntAct ---
    try:
        print("  IntAct 取得中...")
        ia_data = intact.get_interactions(gene_symbol)
        add_edges(ia_data, "IntAct")
        print(f"  IntAct: {len(ia_data)} 件")
    except Exception as e:
        print(f"  [IntAct] エラー: {e}")

    # --- SIGNOR ---
    try:
        print("  SIGNOR 取得中...")
        sg_data = signor.get_interactions(gene_symbol)
        add_edges(sg_data, "SIGNOR")
        print(f"  SIGNOR: {len(sg_data)} 件")
    except Exception as e:
        print(f"  [SIGNOR] エラー: {e}")

    # --- BioGRID (任意) ---
    if use_biogrid:
        try:
            print("  BioGRID 取得中...")
            bg_data = biogrid.get_interactions(gene_symbol, api_key=biogrid_api_key)
            add_edges(bg_data, "BioGRID")
            print(f"  BioGRID: {len(bg_data)} 件")
        except Exception as e:
            print(f"  [BioGRID] エラー: {e}")

    print(f"  ネットワーク: {G.number_of_nodes()} ノード / {G.number_of_edges()} エッジ")
    return G


def _db_color(db: str) -> str:
    return {
        "IntAct":  "#4ECDC4",
        "SIGNOR":  "#45B7D1",
        "BioGRID": "#96CEB4",
    }.get(db, "#DDD")


# ──────────────────────────────────────────────────────────────
# 2. エンリッチメント解析
# ──────────────────────────────────────────────────────────────

def run_network_enrichment(
    G: "nx.Graph",
    top_n: int = 30,
) -> dict:
    """ネットワーク内全遺伝子を g:Profiler でエンリッチメント解析する。

    Returns:
        {
          "gene_list": [...],
          "results":   [enrichment dicts],
          "by_source": {source: [top terms]}
        }
    """
    if G is None:
        return {}

    gene_list = [n for n in G.nodes if n]
    print(f"  エンリッチメント対象: {len(gene_list)} 遺伝子")

    results = enrich_mod.run_enrichment(gene_list, top_n=top_n)
    by_source = enrich_mod.top_terms_by_source(results, top_per_source=5)

    print(f"  エンリッチメント: {len(results)} 有意項目 (FDR<0.05)")
    return {
        "gene_list": gene_list,
        "results":   results,
        "by_source": by_source,
    }


# ──────────────────────────────────────────────────────────────
# 3. 可視化 (pyvis)
# ──────────────────────────────────────────────────────────────

def visualize_network_pyvis(
    G: "nx.Graph",
    gene_symbol: str,
    enrichment: dict | None = None,
    output_path: str = "reports/ppi_network.html",
) -> str | None:
    """pyvis でインタラクティブな PPI ネットワーク HTML を生成する。

    Returns:
        生成した HTML ファイルのパス、または pyvis 未インストール時 None
    """
    try:
        from pyvis.network import Network
    except ImportError:
        print("  [Pyvis] 未インストール: pip install pyvis")
        return None

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # エンリッチメント上位パスウェイでノード色分け
    pathway_gene_map: dict[str, str] = {}
    if enrichment:
        palette = ["#FFD700", "#FF8C00", "#7B68EE", "#20B2AA", "#FF69B4"]
        top_terms = enrichment.get("results", [])[:5]
        for i, term in enumerate(top_terms):
            color = palette[i % len(palette)]
            for g in term.get("genes", []):
                if g.upper() not in pathway_gene_map:
                    pathway_gene_map[g.upper()] = color

    net = Network(
        height="620px", width="100%",
        bgcolor="#1a1a2e", font_color="#eee",
        notebook=True, cdn_resources="in_line",
    )
    net.barnes_hut(spring_length=120, spring_strength=0.05, damping=0.09)

    center = gene_symbol.upper()
    for node, attrs in G.nodes(data=True):
        color = attrs.get("color", "#DDD")
        if node == center:
            color = "#FF6B6B"
        elif node in pathway_gene_map:
            color = pathway_gene_map[node]

        label = node
        db    = attrs.get("db", "")
        title = f"<b>{node}</b><br>DB: {db}"
        size  = 30 if node == center else 16

        net.add_node(node, label=label, color=color, size=size, title=title)

    for src, tgt, attrs in G.edges(data=True):
        weight = attrs.get("weight", 1)
        db     = attrs.get("db", "")
        effect = attrs.get("effect", "")
        net.add_edge(
            src, tgt,
            title=f"DB: {db}<br>Effect: {effect}",
            width=max(1, weight * 1.5),
            color={"color": "#555", "highlight": "#FFF"},
        )

    net.save_graph(output_path)
    return output_path


# ──────────────────────────────────────────────────────────────
# 4. 仮説コンテキスト用サマリー生成
# ──────────────────────────────────────────────────────────────

def network_summary_for_llm(
    G: "nx.Graph",
    gene_symbol: str,
    enrichment: dict,
    max_partners: int = 10,
    max_terms: int = 15,
) -> str:
    """LLM プロンプト用のネットワーク・エンリッチメントサマリーを返す。"""
    if G is None:
        return ""

    center = gene_symbol.upper()
    partners = [n for n in G.neighbors(center)][:max_partners]
    n_nodes  = G.number_of_nodes()
    n_edges  = G.number_of_edges()

    lines = [
        f"## PPI Network ({gene_symbol})",
        f"- Nodes: {n_nodes}, Edges: {n_edges}",
        f"- Key interactors: {', '.join(partners)}",
    ]

    results = (enrichment or {}).get("results", [])[:max_terms]
    if results:
        gprofiler_url = "https://biit.cs.ut.ee/gprofiler/gost"
        lines.append(f"\n## Functional Enrichment (g:Profiler FDR<0.05 — {gprofiler_url})")
        for r in results:
            p = r["p_value"]
            term_id = r.get("term_id", "")
            lines.append(
                f"- [{r['source']}] {r['term_name']} ({term_id}) "
                f"p={p:.2e}, n={r['intersection_size']}"
            )

    return "\n".join(lines)
