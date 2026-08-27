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
import snapshot as snap_mod

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
        r.raise_for_status()
        models = [m["name"] for m in r.json().get("models", [])]
        return jsonify({"available": True, "models": models, "url": OLLAMA_BASE_URL})
    except Exception as e:
        # url とエラー内容を返す。OLLAMA_HOST の設定ミス
        # （例: コンテナ内から localhost を見に行っている等）を
        # ユーザー自身が画面上で気づけるようにするため。
        return jsonify({"available": False, "models": [], "url": OLLAMA_BASE_URL, "error": str(e)})


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
        n_clin = sum(1 for p in result if p.get("is_clinical")) if isinstance(result, list) else 0
        return f"{len(result)} 件の論文 (臨床 {n_clin} 件)"
    if key == "opentargets":
        if isinstance(result, dict):
            score = result.get("association_score") or 0
            dt = result.get("datatype_scores") or {}
            # genetic_association (GWAS/ClinVar等の遺伝的根拠) が無ければ 0。
            # genetic_literature（別のデータタイプ）で代用すると異なる根拠を
            # 「遺伝的」として誤表示してしまうため代用しない。
            gen = dt.get("genetic_association") or 0
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
            if "error" in result:
                return f"エラー: {result['error'][:60]}"
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
            tc = result.get("toxcast") or {}
            if not tc.get("available"):
                return "ToxCast: APIキー未設定"
            return f"ToxCast assays: {tc.get('assay_count', 0)}"
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
                                 "is_clinical": p.get("is_clinical", False),
                                 "pub_types": p.get("pub_types", []),
                                 "abstract": (p.get("abstract", "") or "")[:300]} for p in result]}
        if key == "opentargets" and isinstance(result, dict):
            dt = result.get("datatype_scores") or {}
            drugs = result.get("known_drugs") or []
            assoc_dis = result.get("associated_diseases") or []
            return {
                "ensembl_id": result.get("ensembl_id", ""),
                "disease_id": result.get("disease_id", ""),
                "association_score": result.get("association_score"),
                "datatype_scores": dt,
                # genetic_literature を代用しない（別データタイプの誤表示防止、上の理由と同じ）
                "genetic_score": dt.get("genetic_association"),
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
            return {"accession": result.get("accession", ""),
                    "protein_name": result.get("protein_name", ""),
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
                     "score_inferred": bool(ix.get("score_inferred")),
                     "function": ix.get("partner_function", ""),
                     "protein_name": ix.get("partner_protein_name", ""),
                     "accession": ix.get("partner_accession", "")} for ix in result[:30]]
            return {"interactions": rows}
        if key == "gwas" and isinstance(result, list):
            return {"hits": [{"trait": h.get("trait", ""),
                               "pvalue": h.get("p_value", ""),
                               "variant": ", ".join(h.get("snps", [])),
                               "rsids": h.get("snps", []),
                               "beta": h.get("or_beta", ""),
                               "gwas_url": h.get("gwas_url", "")}
                              for h in result[:10]]}
        if key == "clinvar" and isinstance(result, list):
            return {"total": len(result),
                    "variants": [{"name": v.get("title", "") or v.get("variant_id", ""),
                                   "variant_id": v.get("variant_id", ""),
                                   "significance": v.get("clinical_significance", "") or "—",
                                   "condition": v.get("condition", "") or "—",
                                   "review": v.get("review_status", "")}
                                  for v in result[:100]]}
        if key == "chembl" and isinstance(result, list):
            return {"total": len(result),
                    "drugs": [{"name": d.get("name", "") or d.get("chembl_id", ""),
                                "chembl_id": d.get("chembl_id", ""),
                                "phase": d.get("max_phase", ""),
                                "mechanism": d.get("mechanism", ""),
                                "type": d.get("molecule_type", "")} for d in result[:100]]}
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
                    "key_tissues": result.get("key_tissues", []),
                    "url": result.get("url", "")}
        if key == "hpa" and isinstance(result, dict):
            return {"subcellular": result.get("subcellular", []),
                    "protein_tissue": result.get("protein_tissue", []),
                    "tissues": result.get("tissue_expression", []),
                    "single_cell": result.get("single_cell_expression", []),
                    "cancer": result.get("cancer_expression", []),
                    "url": result.get("url", "")}
        if key == "dgidb" and isinstance(result, list):
            return {"total": len(result),
                    "interactions": [{
                        "drug": d.get("drug_name", ""),
                        # interactionTypes が空の DGIdb レコードが多いため
                        # directionality → sources の順にフォールバックして表示する
                        "type": (d.get("interaction_type") or d.get("directionality") or ""),
                        "approved": d.get("approved", False),
                        "score": d.get("score"),
                        "sources": d.get("sources", []),
                        "pmids": d.get("pmids", []),
                    } for d in result[:100]]}
        if key == "clinicaltrials" and isinstance(result, list):
            return {"total": len(result),
                    "trials": [{"title": t.get("title", "")[:80], "phase": t.get("phase", ""),
                                 "status": t.get("status", ""),
                                 "is_active": t.get("is_active", False),
                                 "nct_id": t.get("nct_id", ""),
                                 "url": t.get("url", ""),
                                 "start_date": t.get("start_date", ""),
                                 # 取得データタブ表示専用（LLMコンテキストには含めない）
                                 "sponsor": t.get("sponsor", ""),
                                 "collaborators": t.get("collaborators", [])} for t in result[:100]]}
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
            tc = result.get("toxcast") or {}
            return {
                "gene": tc.get("gene", ""),
                "toxcast_available": tc.get("available", False),
                "assay_count": tc.get("assay_count", 0),
                "assays": tc.get("assays", []),
                "toxcast_note": tc.get("note", ""),
                "drug_adverse_events": result.get("drug_adverse_events", {}),
            }
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
    string_score = int(ppi.get("string_score", 700))
    min_score    = ppi.get("min_score")
    min_score    = float(min_score) if min_score not in (None, "", "null") else None
    hub_threshold = int(ppi.get("hub_threshold", 1000))
    max_nodes    = int(ppi.get("max_nodes", 100))

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
                "pubmed":         lambda: pubmed.search_pubmed(gene, disease_name, max_results=100, disease_efo_id=disease_id),
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
                "chembl":         lambda: chembl.get_drugs_for_target(gene),
                "gnomad":         lambda: gnomad.get_constraint(gene),
                "gtex":           lambda: gtex.get_tissue_expression(gene),
                "hpa":            lambda: hpa.get_expression_profile(gene),
                "dgidb":          lambda: dgidb.get_interactions(gene),
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

            # clinvar: opentargets の synonyms を使い疾患関連バリアントに絞り込む
            ot_synonyms = ot.get("disease_synonyms", []) if isinstance(ot, dict) else []
            try:
                cv = gwas.get_clinvar_variants(gene, disease_query=disease_name, disease_synonyms=ot_synonyms)
                results["clinvar"] = cv
                send("collector_done", gene=gene, source="clinvar", ok=True,
                     summary=_collector_summary("clinvar", cv, None),
                     data=_collector_data("clinvar", cv))
            except Exception as e:
                errors["clinvar"] = str(e)
                send("collector_done", gene=gene, source="clinvar", ok=False,
                     summary=f"エラー: {e}", data=None)

            # clinicaltrials 再検索（chembl/opentargets の既知薬剤名で intervention 検索を追加）
            # 遺伝子シンボル単独の検索では、薬剤名しか記載されない大半の治験を
            # 取りこぼすため（例: BACE1 阻害薬の治験の多くは "BACE1" と書かれない）
            drug_names = [d.get("name") or d.get("drug") or "" for d in known_drugs]
            try:
                trials = clinicaltrials.get_trials(gene, disease_name, drug_names=drug_names,
                                                    disease_efo_id=disease_id)
                results["clinicaltrials"] = trials
                send("collector_done", gene=gene, source="clinicaltrials", ok=True,
                     summary=_collector_summary("clinicaltrials", trials, None),
                     data=_collector_data("clinicaltrials", trials))
            except Exception as e:
                errors["clinicaltrials"] = str(e)
                send("collector_done", gene=gene, source="clinicaltrials", ok=False,
                     summary=f"エラー: {e}", data=None)

            # pathway_connections + related_gene_papers
            ot_disease_id = ot.get("disease_id") if isinstance(ot, dict) else None
            results["pathway_connections"] = []
            results["related_gene_papers"] = {}
            if ot_disease_id:
                try:
                    send("progress", gene=gene, step="collect",
                         message="パスウェイ隣接遺伝子を解析中...")
                    disease_genes = opentargets.get_disease_top_genes(ot_disease_id, top_n=20)
                    pathway_connections = reactome.find_pathway_connections(
                        gene, disease_genes, max_partners=5)
                    results["pathway_connections"] = pathway_connections

                    # パスウェイ隣接遺伝子の論文を並列取得
                    partners_to_fetch = [
                        c.get("partner", "") for c in pathway_connections[:3]
                        if c.get("partner")
                    ]
                    def _fetch_pp(partner):
                        try:
                            pps = pubmed.search_pubmed(
                                partner, disease_name, max_results=5,
                                disease_efo_id=disease_id)
                            return partner, [p for p in pps if p.get("relevance_score", 0) > 0]
                        except Exception:
                            return partner, []
                    partner_papers: dict = {}
                    if partners_to_fetch:
                        from concurrent.futures import ThreadPoolExecutor as _TPE
                        with _TPE(max_workers=len(partners_to_fetch)) as _ex:
                            for _p, _pp in _ex.map(_fetch_pp, partners_to_fetch):
                                if _pp:
                                    partner_papers[_p] = _pp
                    results["related_gene_papers"] = partner_papers
                except Exception as e:
                    send("progress", gene=gene, step="collect",
                         message=f"パスウェイ解析警告: {e}")

            # Disease pathway enrichment + target fit assessment
            from collectors import gprofiler
            results["pathway_fit"] = {
                "disease_pathways": [], "target_in_disease_pathways": [],
                "pathway_overlap_score": 0.0, "gene_list_size": 0,
            }
            if ot_disease_id:
                try:
                    send("progress", gene=gene, step="collect",
                         message="疾患パスウェイエンリッチメント解析中...")
                    disease_genes_enrich = opentargets.get_disease_top_genes(ot_disease_id, top_n=20)
                    enriched = gprofiler.enrich_gene_list([g["symbol"] for g in disease_genes_enrich])
                    disease_pathway_ids = {p["term_id"] for p in enriched if p["source"] == "REAC"}
                    target_in = reactome.get_gene_pathway_membership(gene, disease_pathway_ids)
                    score = len(target_in) / max(1, min(20, len(disease_pathway_ids)))
                    results["pathway_fit"] = {
                        "disease_pathways": enriched[:20],
                        "target_in_disease_pathways": target_in,
                        "pathway_overlap_score": round(score, 3),
                        "gene_list_size": len(disease_genes_enrich),
                    }
                except Exception as e:
                    send("progress", gene=gene, step="collect",
                         message=f"パスウェイエンリッチメント警告: {e}")

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
                # PPI ネットワーク構築時にスコアが無いエッジへ設定した代用スコア
                # （パートナーの接続数の逆数）を、対応する取得データにも反映する
                inferred_scores = {}
                center_upper = gene.upper()
                for node in ppi_graph.neighbors(center_upper):
                    ed = ppi_graph.edges[center_upper, node]
                    if ed.get("score_inferred"):
                        inferred_scores[node.upper()] = ed.get("score")

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
                        if ix.get("score") is None and p in inferred_scores:
                            ix["score"] = inferred_scores[p]
                            ix["score_inferred"] = True
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
            # コンテキスト長からnum_ctxを動的に決定。
            # 入力 = エビデンス + プロンプトテンプレート（約1500トークン）
            # 出力 = max_tokens（4000）
            # num_ctx < 入力+出力 だとOllamaが無音でプロンプトを打ち切るため
            # 十分な余裕を持たせる。
            ctx_chars = len(context)
            ctx_tokens_est = ctx_chars // 3  # 英語文字÷3でトークン数を粗推定
            template_tokens = 1500            # プロンプトテンプレート分
            generation_tokens = 4000          # max_tokens と合わせる
            num_ctx = max(16384, min(32768, ctx_tokens_est + template_tokens + generation_tokens))
            send("progress", gene=gene, step="llm",
                 message=f"仮説生成中... (コンテキスト約{ctx_tokens_est:,}トークン)")
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
                    num_ctx=num_ctx,
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

            # LLM が誤って自己流の "## References" を書いていれば除去し、
            # サーバー側でエビデンスコンテキストから確定的に生成した
            # リンク付き References に差し替える（Web画面・MD・スナップショット
            # すべてで一貫させるため、ここで hypothesis 本体に結合する）。
            references_section = rpt.references_md(evidence.get("full_references") or {})
            if references_section:
                hypothesis = rpt.strip_llm_references(hypothesis) + "\n\n---\n\n" + references_section
                q.put({"type": "token", "gene": gene, "token": "\n\n---\n\n" + references_section})

            # 使用論文リストをフロントエンドに送信（データ取得タブの Paper タグ対応表示用）
            paper_refs = (evidence.get("full_references") or {}).get("paper", [])
            if paper_refs:
                send("cited_papers", gene=gene, papers=[
                    {"tag": tag, "title": full, "url": url}
                    for tag, full, url in paper_refs
                ])

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
                competitive_section=rpt.competitive_landscape_md(results.get("clinicaltrials")),
                model=model,
            )
            suffix = "JA" if lang == "ja" else "EN"
            rpt_path = pair_dir / f"{ts}_{suffix}.md"
            rpt_path.write_text(md, encoding="utf-8")
            (pair_dir / f"{ts}_raw.json").write_text(
                json.dumps(evidence, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

            # 使用したモデル・PPI設定等のメタ情報（履歴一覧での表示用）
            run_meta = {
                "model": model,
                "lang": lang,
                "disease_id": disease_id,
                "ppi_sources": (["signor"] if use_signor else [])
                               + (["string"] if use_string else [])
                               + (["biogrid"] if use_biogrid_sel else []),
                "string_score": string_score,
                "min_score": min_score,
                "hub_threshold": hub_threshold,
                "max_nodes": max_nodes,
            }
            (pair_dir / f"{ts}_meta.json").write_text(
                json.dumps(run_meta, ensure_ascii=False, indent=2), encoding="utf-8")

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

            # ── 6. Web 表示のスナップショット保存 ──────────────────────────────
            # 取得データ（全ソース）を、Webアプリと同じ見た目でオフライン閲覧
            # できる単一 HTML として書き出す（サーバー不要・ブラウザで直接開ける）。
            collectors_snapshot = {}
            for src in list(COLLECTORS.keys()) + ["toxicity", "clinicaltrials"]:
                result = results.get(src)
                err = errors.get(src)
                collectors_snapshot[src] = {
                    "ok": (err is None and result is not None),
                    "summary": _collector_summary(src, result, err),
                    "data": _collector_data(src, result),
                }
            enr_results_full = (enrichment or {}).get("results", [])[:100]
            excluded_hubs_full = (enrichment or {}).get("excluded_hubs", [])
            snapshot_html = snap_mod.build_snapshot_html(
                gene=gene, disease=disease_name, lang=lang,
                generated_iso=datetime.now().isoformat(),
                collectors=collectors_snapshot,
                hypothesis=hypothesis,
                ppi_image_filename=(Path(ppi_image_rel).name if ppi_image_rel else ""),
                partners=partners,
                partner_functions=partner_functions,
                enrichment_results=enr_results_full,
                excluded_hubs=excluded_hubs_full,
                report_filename=rpt_path.name,
                model=model,
            )
            snapshot_path = pair_dir / f"{ts}_snapshot.html"
            snapshot_path.write_text(snapshot_html, encoding="utf-8")

            send("gene_done", gene=gene, path=str(rpt_path),
                 hypothesis=hypothesis, model=model,
                 ppi_image=ppi_image_rel, partners=partners,
                 partner_functions=partner_functions,
                 enrichment_results=enr_results_full,
                 excluded_hubs=excluded_hubs_full,
                 snapshot_path=str(snapshot_path))

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


# ─── API: history（過去の解析結果の一覧） ──────────────────────────────────
@app.route("/api/history")
def history():
    """reports/ 配下の全遺伝子×疾患ディレクトリをスキャンし、
    保存済みスナップショット/MDレポートの一覧を新しい順で返す。"""
    entries = []
    if REPORTS_DIR.exists():
        for d in REPORTS_DIR.iterdir():
            if not d.is_dir():
                continue
            snapshots = {p.stem.replace("_snapshot", ""): p for p in d.glob("*_snapshot.html")}
            mds       = list(d.glob("*.md"))
            timestamps = set(snapshots.keys())
            for m in mds:
                # ファイル名例: 20260708_220735_JA.md → ts = 20260708_220735
                parts = m.stem.rsplit("_", 1)
                if len(parts) == 2:
                    timestamps.add(parts[0])
            for ts in timestamps:
                snap = snapshots.get(ts)
                md = next((m for m in mds if m.stem.startswith(ts)), None)
                if not snap and not md:
                    continue
                lang = ""
                if md:
                    suffix = md.stem.rsplit("_", 1)[-1]
                    lang = suffix if suffix in ("JA", "EN") else ""
                meta = {}
                meta_path = d / f"{ts}_meta.json"
                if meta_path.exists():
                    try:
                        meta = json.loads(meta_path.read_text(encoding="utf-8"))
                    except Exception:
                        meta = {}
                entries.append({
                    "dir": d.name,
                    "timestamp": ts,
                    "lang": lang,
                    "snapshot_url": f"/reports/{d.name}/{snap.name}" if snap else None,
                    "md_url": f"/reports/{d.name}/{md.name}" if md else None,
                    "model": meta.get("model", ""),
                    "ppi_sources": meta.get("ppi_sources", []),
                    "string_score": meta.get("string_score"),
                    "min_score": meta.get("min_score"),
                    "hub_threshold": meta.get("hub_threshold"),
                    "max_nodes": meta.get("max_nodes"),
                })
    entries.sort(key=lambda e: e["timestamp"], reverse=True)
    return jsonify({"entries": entries[:300]})


@app.route("/")
def index():
    return render_template("index.html")


if __name__ == "__main__":
    port  = int(os.environ.get("PORT", 5000))
    host  = os.environ.get("HOST", "127.0.0.1")
    debug = os.environ.get("FLASK_DEBUG", "1") == "1"
    print("=" * 55)
    print(f"  Drug Hypothesis Generation Web App")
    print(f"  http://{host}:{port}")
    print("=" * 55)
    app.run(debug=debug, threaded=True, host=host, port=port)
