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

def visualize_network_plotly(
    G: "nx.Graph",
    gene_symbol: str,
    enrichment: dict | None = None,
    max_nodes: int = 50,
) -> "go.Figure | None":
    """Plotly でインタラクティブな PPI ネットワーク図を生成する。

    Jupyter でインライン表示可能。pyvis 不要。
    max_nodes: 中心ノード + 上位 N 件（エッジ重み順）に絞る。
    Returns:
        plotly Figure、または依存ライブラリ未インストール時 None
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("  [Plotly] 未インストール: pip install plotly")
        return None

    center = gene_symbol.upper()

    # ── 上位 max_nodes 件に絞ったサブグラフ ─────────────────────────────
    neighbors = sorted(
        G.neighbors(center),
        key=lambda n: G.edges[center, n].get("weight", 1),
        reverse=True,
    )[:max_nodes]
    visible = {center} | set(neighbors)
    subG = G.subgraph(visible)

    n_total = G.number_of_nodes()
    n_shown = subG.number_of_nodes()
    if n_total > n_shown:
        print(f"  表示: {n_shown} ノード（全 {n_total} 件中、上位 {max_nodes} を表示）")

    # ── レイアウト計算 ───────────────────────────────────────────────────
    pos = nx.spring_layout(subG, seed=42, k=2.5 / (n_shown ** 0.5))

    # ── エンリッチメント上位5パスウェイでノード色分け ────────────────────
    PATHWAY_PALETTE = ["#FFD700", "#FF8C00", "#7B68EE", "#20B2AA", "#FF69B4"]
    pathway_gene_map: dict[str, tuple[str, str]] = {}  # node -> (color, term_name)
    if enrichment:
        top_terms = enrichment.get("results", [])[:5]
        for i, term in enumerate(top_terms):
            clr = PATHWAY_PALETTE[i % len(PATHWAY_PALETTE)]
            for g in term.get("genes", []):
                if g.upper() not in pathway_gene_map:
                    pathway_gene_map[g.upper()] = (clr, term["term_name"][:40])

    DB_COLOR = {"IntAct": "#4ECDC4", "SIGNOR": "#45B7D1", "BioGRID": "#96CEB4"}

    # ── エッジトレース ───────────────────────────────────────────────────
    edge_x, edge_y = [], []
    for src, tgt in subG.edges():
        x0, y0 = pos[src]
        x1, y1 = pos[tgt]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        mode="lines",
        line=dict(width=0.8, color="#888"),
        hoverinfo="none",
        showlegend=False,
    )

    # ── ノードトレース（グループ別: center / pathway / db / default） ────
    groups: dict[str, dict] = {
        "center":  {"x": [], "y": [], "text": [], "hover": [], "color": "#FF6B6B", "size": 24, "symbol": "star"},
    }
    for i, (clr, _) in enumerate(dict.fromkeys(pathway_gene_map.values())):
        groups.setdefault(f"pathway_{i}", {"x": [], "y": [], "text": [], "hover": [], "color": clr, "size": 16, "symbol": "circle"})
    for db, clr in DB_COLOR.items():
        groups.setdefault(f"db_{db}", {"x": [], "y": [], "text": [], "hover": [], "color": clr, "size": 12, "symbol": "circle"})
    groups["default"] = {"x": [], "y": [], "text": [], "hover": [], "color": "#AAAAAA", "size": 10, "symbol": "circle"}

    def _assign(node):
        if node == center:
            return "center"
        if node in pathway_gene_map:
            clr, _ = pathway_gene_map[node]
            for i, (c, _) in enumerate(dict.fromkeys(pathway_gene_map.values())):
                if c == clr:
                    return f"pathway_{i}"
        db = subG.nodes[node].get("db", "").split(",")[0]
        if db in DB_COLOR:
            return f"db_{db}"
        return "default"

    for node in subG.nodes():
        x, y = pos[node]
        db = subG.nodes[node].get("db", "")
        degree = subG.degree(node)
        pw_info = pathway_gene_map.get(node, (None, ""))
        hover = (
            f"<b>{node}</b><br>"
            f"DB: {db}<br>"
            f"Degree: {degree}"
            + (f"<br>Pathway: {pw_info[1]}" if pw_info[1] else "")
        )
        grp = _assign(node)
        if grp not in groups:
            grp = "default"
        groups[grp]["x"].append(x)
        groups[grp]["y"].append(y)
        groups[grp]["text"].append(node)
        groups[grp]["hover"].append(hover)

    # ── legend ラベルマッピング ──────────────────────────────────────────
    legend_labels = {"center": f"⬤ {gene_symbol} (target)"}
    if enrichment:
        top_terms = enrichment.get("results", [])[:5]
        seen_colors: dict[str, str] = {}
        for i, (clr, term) in enumerate(dict.fromkeys(pathway_gene_map.values())):
            if clr not in seen_colors:
                seen_colors[clr] = top_terms[i]["term_name"][:35] if i < len(top_terms) else clr
                legend_labels[f"pathway_{i}"] = f"⬤ {seen_colors[clr]}"
    for db in DB_COLOR:
        legend_labels[f"db_{db}"] = f"⬤ {db}"
    legend_labels["default"] = "⬤ other"

    node_traces = []
    for grp, data in groups.items():
        if not data["x"]:
            continue
        node_traces.append(go.Scatter(
            x=data["x"], y=data["y"],
            mode="markers+text",
            marker=dict(
                size=data["size"],
                color=data["color"],
                symbol=data["symbol"],
                line=dict(width=1, color="#222"),
            ),
            text=data["text"],
            textposition="top center",
            textfont=dict(size=9, color="#222"),
            hovertext=data["hover"],
            hoverinfo="text",
            name=legend_labels.get(grp, grp),
        ))

    fig = go.Figure(
        data=[edge_trace] + node_traces,
        layout=go.Layout(
            title=dict(
                text=f"PPI Network — {gene_symbol}  "
                     f"({n_shown} nodes / {subG.number_of_edges()} edges)",
                font=dict(size=15, color="#222"),
                x=0.5,
                xanchor="center",
            ),
            showlegend=True,
            legend=dict(
                font=dict(color="#333", size=10),
                bgcolor="rgba(255,255,255,0.85)",
                bordercolor="#CCC",
                borderwidth=1,
            ),
            hovermode="closest",
            paper_bgcolor="#FFFFFF",
            plot_bgcolor="#FFFFFF",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False,
                       showline=False, mirror=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False,
                       showline=False, mirror=False,
                       scaleanchor="x", scaleratio=1),
            margin=dict(l=20, r=20, t=60, b=20),
            width=650,
            height=650,
        ),
    )
    return fig


# 後方互換エイリアス（旧 pyvis 版が呼ばれた場合でもエラーにならないように）
def visualize_network_pyvis(G, gene_symbol, enrichment=None, output_path="reports/ppi_network.html", max_nodes=30):
    print("⚠ visualize_network_pyvis は非推奨です。visualize_network_plotly を使用してください。")
    fig = visualize_network_plotly(G, gene_symbol, enrichment, max_nodes)
    if fig is not None:
        fig.show()
    return None


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
