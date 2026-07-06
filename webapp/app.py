"""Flask web server for drug hypothesis generation.

Usage:
    cd webapp
    pip install flask
    python app.py

Requires Ollama running: ollama serve
"""
import sys
import json
import time
import queue
import threading
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from flask import Flask, render_template, request, jsonify, Response, send_from_directory
import requests as _requests

import aggregator
import pipeline
import report as rpt
import network as net
import hypothesis as hyp
from llm.ollama_client import OllamaClient
from collectors import opentargets as ot_mod

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
        r = _requests.get("http://localhost:11434/api/tags", timeout=3)
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


# ─── API: analyze (SSE streaming) ─────────────────────────────────────────
@app.route("/api/analyze", methods=["POST"])
def analyze():
    data = request.json or {}
    genes        = data.get("genes", [])
    disease_name = data.get("disease_name", "")
    disease_id   = data.get("disease_id", "")
    lang         = data.get("lang", LANG)
    model        = data.get("model", MODEL)

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
            send("gene_start", gene=gene, index=idx, total=len(genes))

            # 1. collect
            send("progress", gene=gene, step="collecting", message="データ収集中...")
            try:
                evidence = aggregator.collect_all(gene, disease_name, verbose=False, disease_id=disease_id)
            except Exception as e:
                send("gene_error", gene=gene, error=f"データ収集失敗: {e}")
                continue

            # 2-3. PPI + enrichment
            send("progress", gene=gene, step="ppi", message="PPIネットワーク構築中...")
            ppi_graph = None
            enrichment = {}
            try:
                ppi_graph = net.build_ppi_network(gene, use_biogrid=False, use_reactome=True)
                enrichment = net.run_network_enrichment(ppi_graph) if ppi_graph else {}
                if ppi_graph:
                    send("ppi_done", gene=gene,
                         nodes=ppi_graph.number_of_nodes(),
                         edges=ppi_graph.number_of_edges(),
                         partners=net.rank_partners(ppi_graph, gene.upper())[:30])
            except Exception as e:
                send("progress", gene=gene, step="ppi", message=f"PPI警告: {e}")

            # 4. context
            context = aggregator.build_llm_context(evidence, config=None)
            if ppi_graph:
                context += "\n\n" + net.network_summary_for_llm(ppi_graph, gene, enrichment)

            # 5. hypothesis (streaming tokens)
            send("progress", gene=gene, step="llm", message="仮説生成中...")
            hypothesis_parts = []

            def on_token(tok, _gene=gene):
                hypothesis_parts.append(tok)
                q.put({"type": "token", "gene": _gene, "token": tok})

            try:
                hypothesis = hyp.generate_hypothesis(
                    gene, disease_name, context, llm,
                    lang=lang, stream_callback=on_token,
                )
            except Exception as e:
                send("gene_error", gene=gene, error=f"仮説生成失敗: {e}")
                continue

            # 6. save report
            send("progress", gene=gene, step="saving", message="レポート保存中...")
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            pair_dir = REPORTS_DIR / f"{gene}_{disease_name.replace(' ', '_')}"
            pair_dir.mkdir(parents=True, exist_ok=True)

            ppi_section = ""
            ppi_image_rel = ""
            partners = []
            if ppi_graph and ppi_graph.number_of_edges() > 0:
                img_name = f"{ts}_ppi.png"
                partners = net.rank_partners(ppi_graph, gene.upper())[:30]
                if net.render_ppi_image(ppi_graph, gene, str(pair_dir / img_name),
                                        enrichment=enrichment, max_nodes=30):
                    ppi_section = rpt.ppi_md(gene, img_name, partners=partners)
                    ppi_image_rel = str(pair_dir / img_name)

            enr_section = rpt.enrichment_md(enrichment)
            md = rpt.build_report(
                gene, disease_name, lang, hypothesis, context,
                generated_iso=datetime.now().isoformat(),
                ppi_section=ppi_section,
                enrichment_section=enr_section,
            )
            suffix = "JA" if lang == "ja" else "EN"
            rpt_path = pair_dir / f"{ts}_{suffix}.md"
            rpt_path.write_text(md, encoding="utf-8")
            (pair_dir / f"{ts}_raw.json").write_text(
                json.dumps(evidence, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

            send("gene_done", gene=gene, path=str(rpt_path),
                 ppi_image=ppi_image_rel,
                 partners=partners,
                 enrichment_results=(enrichment or {}).get("results", [])[:20])

        send("batch_done", total=len(genes))
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
    import os
    port = int(os.environ.get("PORT", 5000))
    print("=" * 55)
    print(f"  Drug Hypothesis Generation Web App")
    print(f"  http://localhost:{port}")
    print("=" * 55)
    app.run(debug=True, threaded=True, port=port)
