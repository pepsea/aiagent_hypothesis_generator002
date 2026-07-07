"""Flask web server for drug hypothesis generation.

Usage:
    cd webapp
    pip install flask
    python app.py

Requires Ollama running: ollama serve
"""
import os
import sys
import json
import time
import queue
import threading
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def _load_dotenv():
    """プロジェクトルート / webapp の .env を読み込み os.environ に反映する。

    export した環境変数はシェルセッションに紐づき、別プロセス（例:
    サーバーの再起動や別ターミナルから起動した場合）には引き継がれない。
    .env ファイルに書いておけば、起動するたびに確実に読み込まれる。
    既存の環境変数がある場合はそちらを優先し、.env では上書きしない。
    """
    for path in (Path(__file__).parent.parent / ".env", Path(__file__).parent / ".env"):
        if not path.exists():
            continue
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key, value = key.strip(), value.strip().strip('"').strip("'")
            os.environ.setdefault(key, value)


_load_dotenv()

# BIOGRID_API_KEY が設定されていれば BioGRID を PPI に含める（非商用ライセンス）
USE_BIOGRID = bool(os.environ.get("BIOGRID_API_KEY"))

from flask import Flask, render_template, request, jsonify, Response, send_from_directory
import requests as _requests

from concurrent.futures import ThreadPoolExecutor, as_completed

import aggregator
import report as rpt
import network as net
import hypothesis as hyp
from llm.ollama_client import OllamaClient, OLLAMA_BASE_URL
from collectors import (
    pubmed, opentargets, uniprot, intact, gwas, chembl,
    gnomad, gtex, hpa, dgidb, clinicaltrials, alphafold,
    reactome, toxicity, signor, biogrid, string_db,
)

app = Flask(__name__)

# ─── Config ────────────────────────────────────────────────────────────────
MODEL = "qwen2.5:14b"
LANG  = "en"
REPORTS_DIR = Path(__file__).parent.parent / "reports"
REPORTS_DIR.mkdir(exist_ok=True)

_llm = None

def get_llm():
    global _llm
    if _llm is None:
        _llm = OllamaClient(model=MODEL)
    return _llm


# ─── API: config ───────────────────────────────────────────────────────────
@app.route("/api/config", methods=["GET", "POST"])
def config():
    global MODEL, LANG, _llm
    if request.method == "POST":
        data = request.json or {}
        if "model" in data:
            MODEL = data["model"]
            _llm = None  # reset so get_llm() recreates
        if "lang" in data:
            LANG = data["lang"]
    return jsonify({"model": MODEL, "lang": LANG})


@app.route("/api/ollama/status")
def ollama_status():
    try:
        r = _requests.get(f"{OLLAMA_BASE_URL}/api/tags", timeout=3)
        models = [m["name"] for m in r.json().get("models", [])]
        return jsonify({"available": True, "models": models})
    except Exception:
        return jsonify({"available": False, "models": []})


# ─── API: disease search ───────────────────────────────────────────────────
_OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
_DISEASE_SEARCH_Q = """
query($q: String!) {
  search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 15}) {
    hits { id name entity description }
  }
}
"""

@app.route("/api/disease/search", methods=["POST"])
def disease_search():
    q = (request.json or {}).get("query", "").strip()
    if not q:
        return jsonify({"results": []})
    try:
        r = _requests.post(_OT_API, json={"query": _DISEASE_SEARCH_Q, "variables": {"q": q}}, timeout=20)
        hits = [h for h in r.json().get("data", {}).get("search", {}).get("hits", [])
                if h.get("entity") == "disease"]
        return jsonify({"results": [{"id": h["id"], "name": h["name"],
                                      "description": (h.get("description") or "")[:120]} for h in hits]})
    except Exception as e:
        return jsonify({"error": str(e), "results": []})


# ─── API: gene validation ──────────────────────────────────────────────────
@app.route("/api/genes/validate", methods=["POST"])
def validate_genes():
    genes = (request.json or {}).get("genes", [])
    results = []
    for gene in genes:
        gene = gene.strip().upper()
        if not gene:
            continue
        try:
            r = _requests.get(f"https://rest.genenames.org/fetch/symbol/{gene}",
                               headers={"Accept": "application/json"}, timeout=10)
            docs = r.json().get("response", {}).get("docs", [])
            if docs:
                results.append({"symbol": docs[0].get("symbol", gene), "valid": True,
                                 "name": docs[0].get("name", "")})
            else:
                results.append({"symbol": gene, "valid": False, "name": ""})
        except Exception:
            results.append({"symbol": gene, "valid": None, "name": ""})
    return jsonify({"results": results})


# ─── Collector helpers ────────────────────────────────────────────────────
_SRC_LABEL = {
    "pubmed": "PubMed", "opentargets": "OpenTargets", "uniprot": "UniProt",
    "signor": "SIGNOR", "string": "STRING", "biogrid": "BioGRID", "intact": "IntAct",
    "gwas": "GWAS", "clinvar": "ClinVar",
    "chembl": "ChEMBL", "gnomad": "gnomAD", "gtex": "GTEx",
    "hpa": "HPA", "dgidb": "DGIdb", "clinicaltrials": "ClinicalTrials",
    "alphafold": "AlphaFold", "reactome": "Reactome", "toxicity": "Toxicity",
}

def _collector_summary(key: str, result, err: str | None) -> str:
    if err:
        return f"エラー: {err[:80]}"
    if result is None:
        return "データなし"
    if key == "pubmed":
        return f"{len(result)} 件の論文"
    if key == "opentargets":
        if isinstance(result, dict):
            score = result.get("association_score") or 0
            dt = result.get("datatype_scores") or {}
            gen = dt.get("genetic_association") or dt.get("genetic_literature") or 0
            n_drugs = len(result.get("known_drugs") or [])
            return (f"関連スコア {float(score):.3f} / 遺伝的 {float(gen):.3f} / 薬剤 {n_drugs} 件")
        return "取得済み"
    if key == "uniprot":
        if isinstance(result, dict):
            return f"{result.get('protein_name', '')} ({result.get('uniprot_id', '')})"
        return "取得済み"
    if key == "intact":
        if isinstance(result, list):
            partners = list({p for ix in result for p in (ix.get("partners") or [])})
            return f"{len(result)} interactions / {len(partners)} partners"
        return "取得済み"
    if key in ("signor", "biogrid", "string"):
        if isinstance(result, list):
            partners = {ix.get("partner") for ix in result if ix.get("partner")}
            return f"{len(result)} interactions / {len(partners)} partners"
        return "取得済み"
    if key == "gwas":
        return f"{len(result)} ヒット" if isinstance(result, list) else "取得済み"
    if key == "clinvar":
        return f"{len(result)} バリアント" if isinstance(result, list) else "取得済み"
    if key == "chembl":
        return f"{len(result)} 薬剤" if isinstance(result, list) else "取得済み"
    if key == "gnomad":
        if isinstance(result, dict):
            pli = result.get("pLI", result.get("pli"))
            loeuf = result.get("LOEUF", result.get("loeuf"))
            pli_s = f"{pli:.3g}" if isinstance(pli, (int, float)) else "N/A"
            loeuf_s = f"{loeuf:.3g}" if isinstance(loeuf, (int, float)) else "N/A"
            return f"pLI={pli_s}  LOEUF={loeuf_s}"
        return "取得済み"
    if key == "gtex":
        if isinstance(result, dict) and "top_tissues" in result:
            top = result["top_tissues"][0] if result["top_tissues"] else {}
            return f"最高発現: {top.get('tissue', '')} {top.get('tpm', 0):.1f} TPM"
        if isinstance(result, dict) and "error" in result:
            return f"エラー: {result['error'][:60]}"
        return "データなし"
    if key == "hpa":
        if isinstance(result, dict):
            n = len(result.get("tissue_expression", []))
            loc = (result.get("subcellular") or [])
            loc_str = ", ".join(loc[:2]) if loc else "N/A"
            return f"組織数: {n} / 局在: {loc_str}"
        return "取得済み"
    if key == "dgidb":
        return f"{len(result)} 相互作用" if isinstance(result, list) else "取得済み"
    if key == "clinicaltrials":
        return f"{len(result)} 試験" if isinstance(result, list) else "取得済み"
    if key == "alphafold":
        if isinstance(result, dict):
            plddt = result.get("mean_plddt", result.get("plddt"))
            plddt_s = f"{plddt:.1f}" if isinstance(plddt, (int, float)) else "N/A"
            return f"pLDDT={plddt_s}  {result.get('confidence', '')}"
        return "取得済み"
    if key == "reactome":
        return f"{len(result)} パスウェイ" if isinstance(result, list) else "取得済み"
    if key == "toxicity":
        if isinstance(result, dict):
            return f"リスク: {result.get('overall_risk', 'N/A')}"
        return "取得済み"
    return "取得済み"


def _collector_data(key: str, result) -> dict | None:
    """Return a compact, JSON-safe dict for the frontend data tab."""
    if result is None:
        return None
    try:
        if key == "pubmed" and isinstance(result, list):
            return {"papers": [{"title": p.get("title", ""), "journal": p.get("journal", ""),
                                 "year": p.get("year", ""), "pmid": p.get("pmid", ""),
                                 "abstract": (p.get("abstract", "") or "")[:300]} for p in result[:8]]}
        if key == "opentargets" and isinstance(result, dict):
            dt = result.get("datatype_scores") or {}
            drugs = result.get("known_drugs") or []
            assoc_dis = result.get("associated_diseases") or []
            return {
                "association_score": result.get("association_score"),
                "datatype_scores": dt,
                "genetic_score": dt.get("genetic_association") or dt.get("genetic_literature"),
                "drugs": [{"name": d.get("drug","") or d.get("drug_name",""),
                           "phase": d.get("max_phase","") or d.get("maxClinicalStage",""),
                           "indication": d.get("disease","") or d.get("indication","")}
                          for d in drugs[:10]],
                "associated_diseases": [{"name": a.get("disease",""),
                                          "id": a.get("disease_id",""),
                                          "score": a.get("score"),
                                          "datatype_scores": a.get("datatype_scores", {})}
                                        for a in assoc_dis[:20]],
            }
        if key == "uniprot" and isinstance(result, dict):
            # go_terms can be list of str or list of dict
            raw_go = result.get("go_terms", [])[:15]
            go_terms = []
            for g in raw_go:
                if isinstance(g, dict):
                    go_terms.append(g.get("term", g.get("id", "")))
                else:
                    go_terms.append(str(g))
            return {"protein_name": result.get("protein_name", ""),
                    "function": (result.get("function", "") or "")[:500],
                    "subcellular_location": result.get("subcellular_location", [])[:8],
                    "protein_class": result.get("protein_class", []),
                    "go_terms": go_terms}
        if key == "intact" and isinstance(result, list):
            rows = []
            for ix in result[:15]:
                for p in (ix.get("partners") or []):
                    rows.append({"partner": p,
                                 "type": ix.get("interaction_type", ""),
                                 "method": ix.get("detection_method", ""),
                                 "score": ix.get("confidence"),
                                 "pmid": (ix.get("pubmed_ids") or [""])[0]})
            return {"interactions": rows}
        if key in ("signor", "biogrid", "string") and isinstance(result, list):
            rows = [{"partner": ix.get("partner", ""),
                     "effect": ix.get("effect", ""),
                     "mechanism": ix.get("mechanism", ""),
                     "direction": ix.get("direction", ""),
                     "score": ix.get("score"),
                     "function": ix.get("partner_function", ""),
                     "protein_name": ix.get("partner_protein_name", ""),
                     "accession": ix.get("partner_accession", "")} for ix in result[:30]]
            return {"interactions": rows}
        if key == "gwas" and isinstance(result, list):
            return {"hits": [{"trait": h.get("trait", ""),
                               "pvalue": h.get("p_value", ""),
                               "variant": ", ".join(h.get("snps", [])),
                               "beta": h.get("or_beta", "")}
                              for h in result[:10]]}
        if key == "clinvar" and isinstance(result, list):
            return {"variants": [{"name": v.get("title", "") or v.get("variant_id", ""),
                                   "significance": v.get("clinical_significance", "") or "—",
                                   "condition": v.get("condition", "") or "—",
                                   "review": v.get("review_status", "")}
                                  for v in result[:10]]}
        if key == "chembl" and isinstance(result, list):
            return {"drugs": [{"name": d.get("name", "") or d.get("chembl_id", ""),
                                "phase": d.get("max_phase", ""),
                                "mechanism": d.get("mechanism", ""),
                                "type": d.get("molecule_type", "")} for d in result[:10]]}
        if key == "gnomad" and isinstance(result, dict):
            return {"pli": result.get("pLI", result.get("pli")),
                    "loeuf": result.get("LOEUF", result.get("loeuf")),
                    "lof_z": result.get("lof_z"),
                    "obs_lof": result.get("obs_lof"),
                    "exp_lof": result.get("exp_lof"),
                    "essentiality": result.get("essentiality", ""),
                    "url": result.get("url", "")}
        if key == "gtex" and isinstance(result, dict) and "top_tissues" in result:
            return {"tissues": result.get("top_tissues", [])[:10],
                    "key_tissues": result.get("key_tissues", [])}
        if key == "hpa" and isinstance(result, dict):
            tissues = result.get("tissue_expression", [])
            return {"subcellular": result.get("subcellular", []),
                    "protein_tissue": result.get("protein_tissue", []),
                    "tissues": tissues[:15]}
        if key == "dgidb" and isinstance(result, list):
            return {"interactions": [{"drug": d.get("drug_name", ""),
                                       "type": d.get("interaction_type", "")} for d in result[:10]]}
        if key == "clinicaltrials" and isinstance(result, list):
            return {"trials": [{"title": t.get("title", "")[:80], "phase": t.get("phase", ""),
                                 "status": t.get("status", "")} for t in result[:8]]}
        if key == "alphafold" and isinstance(result, dict):
            return {"plddt": result.get("mean_plddt", result.get("plddt")),
                    "confidence": result.get("confidence", ""),
                    "uniprot_id": result.get("uniprot_id", ""),
                    "view_url": result.get("view_url", ""),
                    "pdb_url": result.get("pdb_url", "")}
        if key == "reactome" and isinstance(result, list):
            return {"pathways": [{"name": p.get("name", "") or p.get("pathway_name", ""),
                                   "id": p.get("pathway_id", ""),
                                   "is_disease": p.get("is_disease", False)} for p in result[:15]]}
        if key == "toxicity" and isinstance(result, dict):
            return result
    except Exception:
        pass
    return None


# ─── Job cancellation registry ─────────────────────────────────────────────
_JOBS: dict[str, threading.Event] = {}


class _Cancelled(Exception):
    pass


@app.route("/api/stop", methods=["POST"])
def stop_job():
    job_id = (request.json or {}).get("job_id", "")
    ev = _JOBS.get(job_id)
    if ev:
        ev.set()
        return jsonify({"stopped": True})
    return jsonify({"stopped": False})


# ─── API: analyze (SSE streaming) ─────────────────────────────────────────
@app.route("/api/analyze", methods=["POST"])
def analyze():
    data = request.json or {}
    genes        = data.get("genes", [])
    disease_name = data.get("disease_name", "")
    disease_id   = data.get("disease_id", "")
    lang         = data.get("lang", LANG)
    model        = data.get("model", MODEL)
    job_id       = data.get("job_id", "")
    cancel_ev    = threading.Event()
    if job_id:
        _JOBS[job_id] = cancel_ev

    # ── PPI 設定（ソース選択・クライテリア） ──────────────────────────────
    ppi = data.get("ppi", {}) or {}
    ppi_sources = [s.lower() for s in ppi.get("sources", ["signor"])]
    use_signor  = "signor" in ppi_sources
    use_string  = "string" in ppi_sources
    biogrid_requested = "biogrid" in ppi_sources
    use_biogrid_sel = biogrid_requested and USE_BIOGRID
    string_score = int(ppi.get("string_score", 400))
    min_score    = ppi.get("min_score")
    min_score    = float(min_score) if min_score not in (None, "", "null") else None
    hub_threshold = int(ppi.get("hub_threshold", 3000))
    max_nodes    = int(ppi.get("max_nodes", 30))

    if not genes or not disease_name:
        return jsonify({"error": "genes and disease_name required"}), 400

    # SSE: each event is a JSON line
    q = queue.Queue()

    def send(event_type, **kwargs):
        q.put({"type": event_type, **kwargs})

    def run():
        global _llm
        if model != MODEL or _llm is None:
            _llm = OllamaClient(model=model)
        llm = get_llm()

        send("start", total=len(genes), disease=disease_name)

        for idx, gene in enumerate(genes):
            if cancel_ev.is_set():
                send("stopped", message="ユーザーにより停止されました")
                break
            send("gene_start", gene=gene, index=idx, total=len(genes))

            # ── 1. Parallel data collection with per-collector SSE events ──────
            COLLECTORS = {
                "pubmed":         lambda: pubmed.search_pubmed(gene, disease_name, max_results=8, disease_efo_id=disease_id),
                "opentargets":    lambda: opentargets.get_target_disease_evidence(gene, disease_name, disease_id=disease_id),
                "uniprot":        lambda: uniprot.get_protein_info(gene),
            }
            # 選択された PPI ソースのみ取得データに含める
            if use_signor:
                COLLECTORS["signor"] = lambda: signor.get_interactions(gene)
            if use_string:
                COLLECTORS["string"] = lambda: string_db.get_interactions(gene, required_score=string_score)
            if use_biogrid_sel:
                COLLECTORS["biogrid"] = lambda: biogrid.get_interactions(gene)
            elif biogrid_requested:
                # ユーザーは選択したが BIOGRID_API_KEY 未設定 → 理由を明示して失敗させる
                def _biogrid_missing_key():
                    raise RuntimeError(
                        "BIOGRID_API_KEY が未設定です。サーバー環境変数に設定してください "
                        "(https://webservice.thebiogrid.org/ で無料登録)")
                COLLECTORS["biogrid"] = _biogrid_missing_key
            COLLECTORS.update({
                "gwas":           lambda: gwas.get_gwas_associations(gene, disease_name),
                "clinvar":        lambda: gwas.get_clinvar_variants(gene),
                "chembl":         lambda: chembl.get_drugs_for_target(gene),
                "gnomad":         lambda: gnomad.get_constraint(gene),
                "gtex":           lambda: gtex.get_tissue_expression(gene),
                "hpa":            lambda: hpa.get_expression_profile(gene),
                "dgidb":          lambda: dgidb.get_interactions(gene),
                "clinicaltrials": lambda: clinicaltrials.get_trials(gene, disease_name),
                "alphafold":      lambda: alphafold.get_structure_info(gene),
                "reactome":       lambda: reactome.get_pathways(gene),
            })

            send("collecting_start", gene=gene, sources=list(COLLECTORS.keys()))

            results, errors = {}, {}

            def _run(key, fn):
                try:
                    return key, fn(), None
                except Exception as e:
                    return key, None, str(e)

            with ThreadPoolExecutor(max_workers=10) as ex:
                futures = {ex.submit(_run, k, fn): k for k, fn in COLLECTORS.items()}
                for future in as_completed(futures):
                    key, result, err = future.result()
                    results[key] = result
                    if err:
                        errors[key] = err
                    # summarise for SSE (keep payload small)
                    send("collector_done", gene=gene, source=key,
                         ok=(err is None and result is not None),
                         summary=_collector_summary(key, result, err),
                         data=_collector_data(key, result))

            # toxicity (depends on chembl/opentargets)
            known_drugs = list(results.get("chembl") or [])
            ot = results.get("opentargets")
            if isinstance(ot, dict):
                known_drugs.extend(ot.get("known_drugs", []))
            try:
                tox = toxicity.assess_target_safety(gene, known_drugs)
                results["toxicity"] = tox
                send("collector_done", gene=gene, source="toxicity", ok=True,
                     summary=_collector_summary("toxicity", tox, None),
                     data=_collector_data("toxicity", tox))
            except Exception as e:
                errors["toxicity"] = str(e)
                send("collector_done", gene=gene, source="toxicity", ok=False,
                     summary=f"エラー: {e}", data=None)

            evidence = {"gene": gene, "disease": disease_name,
                        "evidence": results, "collection_errors": errors}

            # ── 2. PPI + enrichment ────────────────────────────────────────────
            send("progress", gene=gene, step="ppi", message="PPIネットワーク構築中...")
            ppi_graph, enrichment = None, {}
            try:
                ppi_graph = net.build_ppi_network(
                    gene,
                    use_signor=use_signor, use_string=use_string,
                    use_biogrid=use_biogrid_sel,
                    string_required_score=string_score, min_score=min_score,
                    use_reactome=False, use_intact=False)
                enrichment = net.run_network_enrichment(
                    ppi_graph, hub_threshold=hub_threshold) if ppi_graph else {}
                if ppi_graph:
                    send("ppi_done", gene=gene,
                         nodes=ppi_graph.number_of_nodes(),
                         edges=ppi_graph.number_of_edges(),
                         partners=net.rank_partners(ppi_graph, gene.upper(), hub_threshold=hub_threshold)[:max_nodes])
            except Exception as e:
                send("progress", gene=gene, step="ppi", message=f"PPI警告: {e}")

            # ── 3. UniProt 機能情報（対象遺伝子 + PPI パートナー） ───────────────
            # 表示用（取得データタブ）: 化合物等は除くがハブは含む「全パートナー」
            #   → タンパク質である限り、ハブでも機能・リンクを表示する
            # 解析用（LLMコンテキスト・エンリッチメント・MDレポート）: ハブも除外
            #   → 仮説生成・機能解析には非特異的なハブを混入させない
            all_functions = {}
            if ppi_graph:
                display_partners = net.rank_partners(
                    ppi_graph, gene.upper(), exclude_hubs=False, exclude_non_gene=True)
                ppi_partners = net.rank_partners(
                    ppi_graph, gene.upper(), hub_threshold=hub_threshold)[:max_nodes]
                send("progress", gene=gene, step="ppi",
                     message="対象遺伝子・PPI遺伝子のUniProt機能情報を取得中...")
                try:
                    all_functions = uniprot.get_functions_for_genes([gene] + display_partners)
                except Exception:
                    all_functions = {}
                # 対象遺伝子の機能は uniprot コレクター結果からも補完
                u = results.get("uniprot")
                if isinstance(u, dict) and gene.upper() not in all_functions:
                    all_functions[gene.upper()] = {
                        "protein_name": u.get("protein_name", ""),
                        "function": u.get("function", ""),
                    }

                # PPI ソース（SIGNOR/STRING/BioGRID）の取得データ表示に
                # 各パートナーの UniProt 機能を付与し、collector_done を再送信して
                # 取得データタブの表示を更新する
                for src in ("signor", "string", "biogrid"):
                    raw = results.get(src)
                    if not isinstance(raw, list) or not raw:
                        continue
                    for ix in raw:
                        p = (ix.get("partner") or "").upper()
                        info = all_functions.get(p)
                        if info:
                            ix["partner_function"] = info.get("function", "")
                            ix["partner_protein_name"] = info.get("protein_name", "")
                            ix["partner_accession"] = info.get("accession", "")
                    send("collector_done", gene=gene, source=src, ok=True,
                         summary=_collector_summary(src, raw, None),
                         data=_collector_data(src, raw))

            # ── 4. LLM context ────────────────────────────────────────────────
            context = aggregator.build_llm_context(evidence, config=None)
            if ppi_graph:
                context += "\n\n" + net.network_summary_for_llm(
                    ppi_graph, gene, enrichment,
                    partner_functions=all_functions, hub_threshold=hub_threshold)

            # ── 4. Hypothesis streaming ───────────────────────────────────────
            send("progress", gene=gene, step="llm", message="仮説生成中...")
            hypothesis_parts = []

            def on_token(tok, _gene=gene):
                if cancel_ev.is_set():
                    raise _Cancelled()   # Ollama ストリーミングを中断
                hypothesis_parts.append(tok)
                q.put({"type": "token", "gene": _gene, "token": tok})

            try:
                hypothesis = hyp.generate_hypothesis(
                    gene, disease_name, context, llm,
                    lang=lang, stream_callback=on_token,
                )
            except _Cancelled:
                send("stopped", message="ユーザーにより停止されました")
                break
            except Exception as e:
                if cancel_ev.is_set():
                    send("stopped", message="ユーザーにより停止されました")
                    break
                send("gene_error", gene=gene, error=f"仮説生成失敗: {e}")
                continue

            # ── 5. Save report ────────────────────────────────────────────────
            send("progress", gene=gene, step="saving", message="レポート保存中...")
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            pair_dir = REPORTS_DIR / f"{gene}_{disease_name.replace(' ', '_')}"
            pair_dir.mkdir(parents=True, exist_ok=True)

            ppi_section, ppi_image_rel, partners = "", "", []
            if ppi_graph and ppi_graph.number_of_edges() > 0:
                img_name = f"{ts}_ppi.png"
                partners = net.rank_partners(ppi_graph, gene.upper(), hub_threshold=hub_threshold)[:max_nodes]
                if net.render_ppi_image(ppi_graph, gene, str(pair_dir / img_name),
                                        enrichment=enrichment, max_nodes=max_nodes):
                    ppi_section = rpt.ppi_md(gene, img_name, partners=partners,
                                             functions=all_functions)
                    ppi_image_rel = str(pair_dir / img_name)

            md = rpt.build_report(
                gene, disease_name, lang, hypothesis, context,
                generated_iso=datetime.now().isoformat(),
                ppi_section=ppi_section,
                enrichment_section=rpt.enrichment_md(enrichment),
            )
            suffix = "JA" if lang == "ja" else "EN"
            rpt_path = pair_dir / f"{ts}_{suffix}.md"
            rpt_path.write_text(md, encoding="utf-8")
            (pair_dir / f"{ts}_raw.json").write_text(
                json.dumps(evidence, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

            # 機能情報を {gene, protein_name, function} のリストに整形（対象を先頭）
            func_order = [gene.upper()] + [p for p in partners if p.upper() != gene.upper()]
            partner_functions = []
            _seen = set()
            for g in func_order:
                info = all_functions.get(g.upper())
                if not info or g.upper() in _seen or not info.get("function"):
                    continue
                _seen.add(g.upper())
                partner_functions.append({
                    "gene": g,
                    "protein_name": info.get("protein_name", ""),
                    "function": info.get("function", ""),
                    "is_target": g.upper() == gene.upper(),
                })

            send("gene_done", gene=gene, path=str(rpt_path),
                 hypothesis=hypothesis,
                 ppi_image=ppi_image_rel, partners=partners,
                 partner_functions=partner_functions,
                 enrichment_results=(enrichment or {}).get("results", [])[:100],
                 excluded_hubs=(enrichment or {}).get("excluded_hubs", []))

        if not cancel_ev.is_set():
            send("batch_done", total=len(genes))
        _JOBS.pop(job_id, None)
        q.put(None)  # sentinel

    thread = threading.Thread(target=run, daemon=True)
    thread.start()

    def event_stream():
        while True:
            item = q.get()
            if item is None:
                break
            yield f"data: {json.dumps(item, ensure_ascii=False)}\n\n"

    return Response(event_stream(), mimetype="text/event-stream",
                    headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})


# ─── Serve report files (images, markdown) ────────────────────────────────
@app.route("/reports/<path:filename>")
def serve_report(filename):
    return send_from_directory(str(REPORTS_DIR), filename)


@app.route("/")
def index():
    return render_template("index.html")


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    print("=" * 55)
    print(f"  Drug Hypothesis Generation Web App")
    print(f"  http://localhost:{port}")
    print("=" * 55)
    app.run(debug=True, threaded=True, port=port)
