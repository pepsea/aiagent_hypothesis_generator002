"""PPI ネットワーク構築 + エンリッチメント解析統合モジュール.

データソース: IntAct (EMBL-EBI) + SIGNOR (UNIROMA2) + Reactome（すべて CC BY 4.0・商用可）
             BioGRID は非商用ライセンスのため任意（デフォルト無効）
エンリッチメント: g:Profiler (GO/Reactome/WikiPathways — KEGG等の商用ソースは不使用)
可視化: pyvis (インタラクティブHTML) または networkx/matplotlib
"""
from __future__ import annotations

import os
import json
import time
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from typing import Optional

try:
    import networkx as nx
    HAS_NX = True
except ImportError:
    HAS_NX = False

from collectors import intact, signor, biogrid, string_db, enrichment as enrich_mod
from collectors import reactome as reactome_mod

import re

# ── ハブ遺伝子の判定基準（2つの客観ルールの OR）───────────────────────────
# (1) グローバル相互作用数（IntAct 実測値）が閾値を超える → 非特異的ハブ
# (2) 既知の「無差別に相互作用するタンパク質族」に遺伝子名が一致する
#     （G タンパク質サブユニット、ユビキチン、チューブリン等）。これらは
#     相互作用数が中程度でもシグナル/構造の非特異ハブとして機能するため、
#     数だけでは捕捉できない。族はパターン（正規表現）で定義し再現可能。
HUB_DEGREE_THRESHOLD = int(os.environ.get("HUB_DEGREE_THRESHOLD", "3000"))

# 非特異ハブ族の遺伝子名パターン（HGNC 命名規則に基づく）
HUB_FAMILY_PATTERNS = [
    r"^GNA[0-9A-Z]",  # G タンパク質 α サブユニット (GNAI2, GNAQ, GNAS, GNA11, GNAO1 …)
    r"^GNB[0-9]",     # G タンパク質 β サブユニット (GNB1, GNB2 …)
    r"^GNG[0-9]",     # G タンパク質 γ サブユニット (GNG2 …)
    r"^UB[ABC][0-9]?$",  # ユビキチン (UBC, UBB, UBA52 …)
    r"^RPS27A$",      # ユビキチン融合
    r"^TUB[AB]",      # チューブリン (TUBB, TUBA1A, TUBB4B …)
    r"^ACT[BG][0-9]?$",  # アクチン (ACTB, ACTG1)
    r"^HSP(90|A)",    # 主要シャペロン (HSP90AA1, HSPA8 …)
    r"^YWHA[BEGHQZ]$",  # 14-3-3 (YWHAZ, YWHAE …)
    r"^SUMO[0-9]$",   # SUMO
]
_HUB_FAMILY_RE = re.compile("|".join(HUB_FAMILY_PATTERNS))


def is_hub_family(gene_symbol: str) -> bool:
    """既知の非特異ハブ族（G タンパク質・ユビキチン等）に一致するか。"""
    return bool(_HUB_FAMILY_RE.match(gene_symbol.upper()))

_HUB_CACHE_DIR = Path(__file__).parent / "ppi_cache" / "hub_degree"
_HUB_CACHE_TTL = 30 * 24 * 3600   # 30日（相互作用数は頻繁には変わらない）
_hub_mem_cache: dict[str, int] = {}


def global_interactor_count(gene_symbol: str) -> int:
    """IntAct における遺伝子のグローバル相互作用数を返す（客観的なハブ指標）。

    取得失敗時は -1 を返す（除外判定に使わない）。結果はメモリ+ファイルにキャッシュ。
    """
    key = gene_symbol.upper()
    if key in _hub_mem_cache:
        return _hub_mem_cache[key]

    cache_file = _HUB_CACHE_DIR / f"{key}.json"
    if cache_file.exists() and time.time() - cache_file.stat().st_mtime < _HUB_CACHE_TTL:
        try:
            n = json.loads(cache_file.read_text())["count"]
            _hub_mem_cache[key] = n
            return n
        except Exception:
            pass

    try:
        r = requests.get(
            f"https://www.ebi.ac.uk/intact/ws/interaction/findInteractions/{gene_symbol}",
            params={"page": 0, "pageSize": 1, "query": "species:9606"}, timeout=15)
        r.raise_for_status()
        n = int(r.json().get("totalElements", -1))
    except Exception:
        n = -1

    if n >= 0:
        _HUB_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps({"count": n}), encoding="utf-8")
    _hub_mem_cache[key] = n
    return n


# ──────────────────────────────────────────────────────────────
# 1. ネットワーク構築
# ──────────────────────────────────────────────────────────────

def build_ppi_network(
    gene_symbol: str,
    use_biogrid: bool = True,
    biogrid_api_key: str | None = None,
    use_reactome: bool = False,
    use_intact: bool = False,
    use_signor: bool = True,
    use_string: bool = False,
    string_required_score: int = 400,
    min_score: float | None = None,
) -> Optional["nx.Graph"]:
    """選択された PPI ソースから NetworkX グラフを構築する。

    ソース選択: use_signor / use_string / use_biogrid /（use_intact / use_reactome）。
    string_required_score: STRING の信頼度閾値（0–1000、400=中, 700=高）。
    min_score: 全ソース共通のエッジスコア下限（None なら適用しない）。
    BioGRID は非商用・学術利用限定ライセンス（BIOGRID_API_KEY が必要）。

    Returns:
        nx.Graph: ノード属性に db, color, entity_type を持つグラフ
                  networkx 未インストールの場合は None
    """
    if not HAS_NX:
        print("  [Network] networkx が未インストールです: pip install networkx")
        return None

    G = nx.Graph()
    center = gene_symbol.upper()
    G.add_node(center, color="#FF6B6B", size=25, db="center", direct_partner=True)

    def _item_partners(item: dict):
        """1件の相互作用から (partner, score) のリストを返す。

        2形式に対応:
          - {"source": ..., "target": ..., "score"/"confidence": ...}  (SIGNOR/Reactome)
          - {"partners": [...], "confidence": ...}                      (IntAct)
        """
        score = item.get("score", item.get("confidence", None))
        try:
            score = float(score) if score is not None else None
        except (TypeError, ValueError):
            score = None

        ptype = (item.get("partner_type") or "gene").strip().lower()
        out = []
        src = (item.get("source") or "").strip().upper()
        tgt = (item.get("target") or "").strip().upper()
        if src and tgt:
            partner = tgt if src == center else src
            if partner and partner != center:
                out.append((partner, score, ptype))
        for p in item.get("partners", []) or []:
            p = str(p).strip().upper()
            if p and p != center:
                out.append((p, score, ptype))
        return out

    def add_edges(interactions: list[dict], source_label: str):
        for item in interactions:
            for partner, score, ptype in _item_partners(item):
                # ノード追加
                if partner not in G:
                    G.add_node(partner, color=_db_color(source_label),
                               size=15, db=source_label, direct_partner=True,
                               entity_type=ptype)
                else:
                    if source_label not in G.nodes[partner].get("db", ""):
                        G.nodes[partner]["db"] += f",{source_label}"
                    # いずれかのソースで gene と判定されれば gene を優先
                    if ptype == "gene":
                        G.nodes[partner]["entity_type"] = "gene"

                # エッジ追加 or 属性更新
                if G.has_edge(center, partner):
                    ed = G.edges[center, partner]
                    ed["weight"] = ed.get("weight", 1) + 1
                    ed.setdefault("dbs", set()).add(source_label)
                    if source_label not in ed.get("db", ""):
                        ed["db"] += f",{source_label}"
                    if score is not None:
                        prev = ed.get("score")
                        ed["score"] = score if prev is None else max(prev, score)
                else:
                    G.add_edge(center, partner,
                               weight=1,
                               effect=item.get("effect", ""),
                               mechanism=item.get("mechanism", ""),
                               db=source_label,
                               dbs={source_label},
                               score=score)

    # --- IntAct（任意） ---
    if use_intact:
        try:
            print("  IntAct 取得中...")
            ia_data = intact.get_interactions(gene_symbol)
            add_edges(ia_data, "IntAct")
            print(f"  IntAct: {len(ia_data)} 件")
        except Exception as e:
            print(f"  [IntAct] エラー: {e}")

    # --- SIGNOR ---
    if use_signor:
        try:
            print("  SIGNOR 取得中...")
            sg_data = signor.get_interactions(gene_symbol)
            add_edges(sg_data, "SIGNOR")
            print(f"  SIGNOR: {len(sg_data)} 件")
        except Exception as e:
            print(f"  [SIGNOR] エラー: {e}")

    # --- STRING（任意、信頼度閾値付き） ---
    if use_string:
        try:
            print(f"  STRING 取得中... (required_score={string_required_score})")
            st_data = string_db.get_interactions(gene_symbol, required_score=string_required_score)
            add_edges(st_data, "STRING")
            print(f"  STRING: {len(st_data)} 件")
        except Exception as e:
            print(f"  [STRING] エラー: {e}")

    # --- BioGRID (任意) ---
    if use_biogrid:
        try:
            print("  BioGRID 取得中...")
            bg_data = biogrid.get_interactions(gene_symbol, api_key=biogrid_api_key)
            add_edges(bg_data, "BioGRID")
            print(f"  BioGRID: {len(bg_data)} 件")
        except Exception as e:
            print(f"  [BioGRID] エラー: {e}")

    # --- Reactome ---
    if use_reactome:
        try:
            print("  Reactome 取得中...")
            rc_data = reactome_mod.get_interactions(gene_symbol)
            add_edges(rc_data, "Reactome")
            print(f"  Reactome: {len(rc_data)} 件")
        except Exception as e:
            print(f"  [Reactome] エラー: {e}")

    # --- スコアが無いエッジへのフォールバック ---
    # PPI データベースによってはエッジにスコア（信頼度）を提供しない場合がある
    # （例: BioGRID の SCORE は多くのレコードで空）。この場合、パートナー遺伝子の
    # グローバル接続数（IntAct 実測の総相互作用数）の逆数を代用スコアとする。
    # 接続数が少ない（無差別なハブでない）パートナーほど、その相互作用が特異的
    # ＝意味がある可能性が高いという考え方に基づく。
    no_score_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("score") is None]
    if no_score_edges:
        partners_needing_count = [v if u == center else u for u, v in no_score_edges]
        with ThreadPoolExecutor(max_workers=10) as ex:
            counts = dict(zip(partners_needing_count,
                              ex.map(global_interactor_count, partners_needing_count)))
        n_inferred = 0
        for u, v in no_score_edges:
            partner = v if u == center else u
            count = counts.get(partner, -1)
            if count and count > 0:
                G.edges[u, v]["score"] = 1.0 / count
                G.edges[u, v]["score_inferred"] = True
                n_inferred += 1
        if n_inferred:
            print(f"  スコア未提供の {n_inferred} エッジに接続数逆数を代用スコアとして設定")

    # --- スコア下限フィルタ（共通クライテリア） ---
    if min_score is not None:
        drop = [(u, v) for u, v, d in G.edges(data=True)
                if d.get("score") is not None and d["score"] < min_score]
        G.remove_edges_from(drop)
        # 孤立ノード（中心以外）を削除
        isolated = [n for n in list(G.nodes)
                    if n != center and G.degree(n) == 0]
        G.remove_nodes_from(isolated)
        if drop:
            print(f"  スコア下限 {min_score} で {len(drop)} エッジ除外")

    print(f"  ネットワーク: {G.number_of_nodes()} ノード / {G.number_of_edges()} エッジ")
    return G


def _db_color(db: str) -> str:
    return {
        "IntAct":   "#4ECDC4",
        "SIGNOR":   "#45B7D1",
        "STRING":   "#B39DDB",
        "Reactome": "#FFB347",
        "BioGRID":  "#96CEB4",
    }.get(db, "#DDD")


def _n_distinct_dbs(edge: dict) -> int:
    """エッジが何種類のDBで裏付けられているか。"""
    dbs = edge.get("dbs")
    if not dbs:
        dbs = {d for d in (edge.get("db", "") or "").split(",") if d}
    return len(dbs)


def rank_partners(
    G: "nx.Graph",
    center: str,
    exclude_hubs: bool = True,
    exclude_non_gene: bool = True,
    hub_threshold: int | None = None,
) -> list[str]:
    """PPIパートナーを優先度順に並べて返す。

    優先順位（降順）:
      1. 複数DBで共通する遺伝子（裏付けDB数が多いほど上位）
      2. スコアが高いもの（IntAct intactScore / SIGNOR score など）
      3. エッジの重み（観測された相互作用の回数）

    exclude_hubs: True の場合、ハブ遺伝子（既知の非特異ハブ族、または
        グローバル相互作用数が hub_threshold 超）をリストから除外する。
        仮説生成用コンテキスト（LLM・エンリッチメント・MDレポート）に
        非特異的なハブが紛れ込むのを防ぐ。エンリッチメント解析の
        除外基準（run_network_enrichment）と同じ判定を用いる。
        取得データタブでの機能表示など「全パートナーを見せたいが解析対象
        からは外したい」場合は False にする（exclude_non_gene は独立）。
    exclude_non_gene: True の場合、化合物・フェノタイプ等の非遺伝子ノードを
        常に除外する（exclude_hubs の設定に関わらず適用）。
    """
    def key(n):
        ed = G.edges[center, n]
        n_db  = _n_distinct_dbs(ed)
        score = ed.get("score")
        score = score if score is not None else -1.0
        weight = ed.get("weight", 1)
        return (n_db, score, weight)

    ranked = sorted(G.neighbors(center), key=key, reverse=True)
    candidates = [n for n in ranked if n != center]

    excluded = set()
    if exclude_non_gene:
        excluded |= {n for n in candidates
                     if G.nodes[n].get("entity_type", "gene") != "gene"}

    if exclude_hubs:
        threshold = hub_threshold if hub_threshold is not None else HUB_DEGREE_THRESHOLD
        family_hubs = {n for n in candidates if is_hub_family(n)}
        # 残った候補（族に該当しない・既に除外対象でない）だけ次数を問い合わせる
        to_check = [n for n in candidates if n not in family_hubs and n not in excluded]
        degree_hubs = set()
        if to_check:
            with ThreadPoolExecutor(max_workers=10) as ex:
                counts = dict(zip(to_check, ex.map(global_interactor_count, to_check)))
            degree_hubs = {n for n, c in counts.items() if c > threshold}
        excluded |= family_hubs | degree_hubs

    return [n for n in ranked if n not in excluded]


# ──────────────────────────────────────────────────────────────
# 2. エンリッチメント解析
# ──────────────────────────────────────────────────────────────

def run_network_enrichment(
    G: "nx.Graph",
    top_n: int = 30,
    exclude_hubs: bool = True,
    hub_threshold: int | None = None,
) -> dict:
    """ネットワーク内全遺伝子を g:Profiler でエンリッチメント解析する。

    exclude_hubs: True の場合、ハブ遺伝子を機能解析から除外。
        ハブ判定は IntAct のグローバル相互作用数が hub_threshold を超えるか、
        で客観的に行う（恣意的なリストではなく実測値ベース）。
    hub_threshold: 相互作用数の閾値（None なら HUB_DEGREE_THRESHOLD）。

    Returns:
        {
          "gene_list":     [...],
          "results":       [enrichment dicts],
          "by_source":     {source: [top terms]},
          "excluded_hubs": [{"gene", "interactor_count"}],
        }
    """
    if G is None:
        return {}

    threshold = hub_threshold if hub_threshold is not None else HUB_DEGREE_THRESHOLD

    # 遺伝子/タンパク質ノードのみを対象にする（化合物・phenotype 等を除外）
    # entity_type が無い中心ノード等は gene とみなす
    gene_list = [n for n in G.nodes
                 if n and G.nodes[n].get("entity_type", "gene") == "gene"]
    excluded = [n for n in G.nodes if n and G.nodes[n].get("entity_type", "gene") != "gene"]
    if excluded:
        print(f"  エンリッチメント除外（非遺伝子）: {', '.join(excluded)}")

    # ハブ遺伝子を除外（中心遺伝子は必ず残す）
    #   ルール(1): グローバル相互作用数 > threshold
    #   ルール(2): 既知の非特異ハブ族（G タンパク質・ユビキチン等）に一致
    excluded_hubs = []
    if exclude_hubs and gene_list:
        center = next((n for n in G.nodes if G.nodes[n].get("db") == "center"), None)
        candidates = [g for g in gene_list if g != center]
        with ThreadPoolExecutor(max_workers=10) as ex:
            counts = dict(zip(candidates, ex.map(global_interactor_count, candidates)))
        hubs = {}
        for g in candidates:
            if is_hub_family(g):
                hubs[g] = {"count": counts.get(g, -1), "reason": "family"}
            elif counts.get(g, -1) > threshold:
                hubs[g] = {"count": counts[g], "reason": "degree"}
        if hubs:
            order = sorted(hubs, key=lambda x: hubs[x]["count"], reverse=True)
            excluded_hubs = [{"gene": g, "interactor_count": hubs[g]["count"],
                              "reason": hubs[g]["reason"]} for g in order]
            desc = ", ".join(
                f"{g}({hubs[g]['count'] if hubs[g]['count'] >= 0 else '?'}"
                f"{'/族' if hubs[g]['reason']=='family' else ''})" for g in order)
            print(f"  エンリッチメント除外（ハブ）: {desc}")
            gene_list = [g for g in gene_list if g not in hubs]

    print(f"  エンリッチメント対象: {len(gene_list)} 遺伝子")

    results = enrich_mod.run_enrichment(gene_list, top_n=top_n)
    by_source = enrich_mod.top_terms_by_source(results, top_per_source=5)

    print(f"  エンリッチメント: {len(results)} 有意項目 (FDR<0.05)")
    return {
        "gene_list":     gene_list,
        "results":       results,
        "by_source":     by_source,
        "excluded_hubs": excluded_hubs,
    }


# ──────────────────────────────────────────────────────────────
# 3. 可視化（静的 PNG）
# ──────────────────────────────────────────────────────────────

def render_ppi_image(
    G: "nx.Graph",
    gene_symbol: str,
    out_path: str,
    enrichment: dict | None = None,
    max_nodes: int = 30,
    dpi: int = 130,
) -> str | None:
    """PPI ネットワークを静的 PNG として保存する（レポート埋め込み用）。

    レイアウト: 対象遺伝子を上部中央、PPI パートナーを下部に散らして配置。
    Returns: 保存した PNG パス（成功時）、None（依存欠如・データなし時）。
    """
    if G is None or not HAS_NX:
        return None
    try:
        import math
        import matplotlib
        matplotlib.use("Agg")  # GUI 不要のバックエンド
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [PPI image] matplotlib 未インストール: pip install matplotlib")
        return None

    center = gene_symbol.upper()
    if center not in G:
        return None

    # ── 上位パートナーを抽出（DB共通数→スコア→重み順、最大 max_nodes 件） ──
    neighbors = rank_partners(G, center)[:max_nodes]
    if not neighbors:
        return None

    # ── レイアウト: 中心を上部、パートナーを下部に散らす ─────────────
    pos = {center: (0.5, 1.06)}
    n = len(neighbors)
    cols = max(1, math.ceil(math.sqrt(n * 1.8)))   # 横長に散らす
    rows = math.ceil(n / cols)
    for i, node in enumerate(neighbors):
        r, c = divmod(i, cols)
        # 列方向に均等配置 + 行ごとに段差、わずかなジッターでばらけさせる
        x = (c + 0.5) / cols
        jitter = 0.04 * (1 if (i % 2 == 0) else -1)
        y = 0.55 - (r / max(1, rows)) * 0.55 + jitter
        pos[node] = (x, y)

    # ── DB別の色分け（複数DB共通はハイライト） ──────────────────────
    MULTI_COLOR = "#8E44AD"  # 2つ以上のDBで共通 = 高信頼

    def node_db_color(node: str) -> str:
        ed = G.edges[center, node]
        dbs = ed.get("dbs") or {d for d in (ed.get("db", "") or "").split(",") if d}
        if len(dbs) >= 2:
            return MULTI_COLOR
        only = next(iter(dbs)) if dbs else ""
        return _db_color(only)

    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.15, 1.20)
    ax.axis("off")

    # エッジ（中心 → 各パートナー）
    cx, cy = pos[center]
    for node in neighbors:
        nx_, ny_ = pos[node]
        ax.plot([cx, nx_], [cy, ny_], color="#C8C8C8", lw=0.8, zorder=1)

    # パートナーノード
    used_dbs = set()
    has_multi = False
    for node in neighbors:
        x, y = pos[node]
        ed = G.edges[center, node]
        dbs = ed.get("dbs") or {d for d in (ed.get("db", "") or "").split(",") if d}
        if len(dbs) >= 2:
            has_multi = True
        else:
            used_dbs |= dbs
        ax.scatter([x], [y], s=320, c=node_db_color(node), edgecolors="#222",
                   linewidths=0.8, zorder=2)
        ax.text(x, y - 0.055, node, ha="center", va="top",
                fontsize=7.5, color="#222", zorder=3)

    # 中心ノード（星）
    ax.scatter([cx], [cy], s=900, c="#FF6B6B", marker="*",
               edgecolors="#222", linewidths=1.0, zorder=4)
    ax.text(cx, cy + 0.06, center, ha="center", va="bottom",
            fontsize=12, fontweight="bold", color="#B22222", zorder=5)

    # 凡例（DBソース）
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="*", color="w", label=f"{gene_symbol} (target)",
                      markerfacecolor="#FF6B6B", markeredgecolor="#222", markersize=15)]
    for db in ("IntAct", "SIGNOR", "Reactome", "BioGRID"):
        if db in used_dbs:
            handles.append(Line2D([0], [0], marker="o", color="w", label=db,
                                  markerfacecolor=_db_color(db), markeredgecolor="#222", markersize=9))
    if has_multi:
        handles.append(Line2D([0], [0], marker="o", color="w", label="multiple DBs",
                              markerfacecolor=MULTI_COLOR, markeredgecolor="#222", markersize=9))
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True,
              framealpha=0.9, edgecolor="#CCC")

    ax.set_title(
        f"PPI Network — {gene_symbol}  "
        f"({G.number_of_nodes()} nodes / {G.number_of_edges()} edges; "
        f"showing top {len(neighbors)})",
        fontsize=11, color="#222",
    )

    import os
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


# ──────────────────────────────────────────────────────────────
# 4. 仮説コンテキスト用サマリー生成
# ──────────────────────────────────────────────────────────────

def network_summary_for_llm(
    G: "nx.Graph",
    gene_symbol: str,
    enrichment: dict,
    max_partners: int = 10,
    max_terms: int = 15,
    partner_functions: dict | None = None,
    hub_threshold: int | None = None,
) -> str:
    """LLM プロンプト用のネットワーク・エンリッチメントサマリーを返す。

    partner_functions: {GENE(upper): {"protein_name","function"}} を渡すと
        各 PPI パートナーの UniProt 機能情報を仮説生成コンテキストに含める。
    ハブ遺伝子は rank_partners の既定どおり除外される（仮説生成対象から外す）。
    """
    if G is None:
        return ""

    center = gene_symbol.upper()
    partners = rank_partners(G, center, hub_threshold=hub_threshold)[:max_partners]
    n_nodes  = G.number_of_nodes()
    n_edges  = G.number_of_edges()

    lines = [
        f"## PPI Network ({gene_symbol})",
        f"- Nodes: {n_nodes}, Edges: {n_edges}",
        f"- Key interactors: {', '.join(partners)}",
    ]

    # PPI パートナーの UniProt 機能情報
    if partner_functions:
        lines.append("\n## PPI Partner Functions (UniProt)")
        for p in partners:
            info = partner_functions.get(p.upper())
            if not info or not info.get("function"):
                continue
            pname = info.get("protein_name", "")
            func = info["function"]
            if len(func) > 400:
                func = func[:400].rsplit(" ", 1)[0] + " ..."
            label = f"{p} ({pname})" if pname else p
            lines.append(f"- **{label}**: {func}")

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
