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
    reactome, toxicity,
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
    "intact": "IntAct", "gwas": "GWAS", "clinvar": "ClinVar",
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

            # ── 1. Parallel data collection with per-collector SSE events ──────
            COLLECTORS = {
                "pubmed":         lambda: pubmed.search_pubmed(gene, disease_name, max_results=8, disease_efo_id=disease_id),
                "opentargets":    lambda: opentargets.get_target_disease_evidence(gene, disease_name, disease_id=disease_id),
                "uniprot":        lambda: uniprot.get_protein_info(gene),
                "intact":         lambda: intact.get_interactions(gene, max_results=15),
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
            }

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
                ppi_graph = net.build_ppi_network(gene, use_biogrid=USE_BIOGRID, use_reactome=True)
                enrichment = net.run_network_enrichment(ppi_graph) if ppi_graph else {}
                if ppi_graph:
                    send("ppi_done", gene=gene,
                         nodes=ppi_graph.number_of_nodes(),
                         edges=ppi_graph.number_of_edges(),
                         partners=net.rank_partners(ppi_graph, gene.upper())[:30])
            except Exception as e:
                send("progress", gene=gene, step="ppi", message=f"PPI警告: {e}")

            # ── 3. LLM context ────────────────────────────────────────────────
            context = aggregator.build_llm_context(evidence, config=None)
            if ppi_graph:
                partners_for_fn = net.rank_partners(ppi_graph, gene.upper())[:10]
                partner_fns = {}
                if partners_for_fn:
                    send("progress", gene=gene, step="ppi",
                         message="PPIパートナーのUniProt機能情報を取得中...")
                    try:
                        partner_fns = uniprot.get_functions_for_genes(partners_for_fn)
                    except Exception:
                        partner_fns = {}
                context += "\n\n" + net.network_summary_for_llm(
                    ppi_graph, gene, enrichment, partner_functions=partner_fns)

            # ── 4. Hypothesis streaming ───────────────────────────────────────
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

            # ── 5. Save report ────────────────────────────────────────────────
            send("progress", gene=gene, step="saving", message="レポート保存中...")
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            pair_dir = REPORTS_DIR / f"{gene}_{disease_name.replace(' ', '_')}"
            pair_dir.mkdir(parents=True, exist_ok=True)

            ppi_section, ppi_image_rel, partners = "", "", []
            if ppi_graph and ppi_graph.number_of_edges() > 0:
                img_name = f"{ts}_ppi.png"
                partners = net.rank_partners(ppi_graph, gene.upper())[:30]
                if net.render_ppi_image(ppi_graph, gene, str(pair_dir / img_name),
                                        enrichment=enrichment, max_nodes=30):
                    ppi_section = rpt.ppi_md(gene, img_name, partners=partners)
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

            send("gene_done", gene=gene, path=str(rpt_path),
                 hypothesis=hypothesis,
                 ppi_image=ppi_image_rel, partners=partners,
                 enrichment_results=(enrichment or {}).get("results", [])[:100],
                 excluded_hubs=(enrichment or {}).get("excluded_hubs", []))

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
    port = int(os.environ.get("PORT", 5000))
    print("=" * 55)
    print(f"  Drug Hypothesis Generation Web App")
    print(f"  http://localhost:{port}")
    print("=" * 55)
    app.run(debug=True, threaded=True, port=port)
