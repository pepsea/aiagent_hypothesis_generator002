"""Aggregate multi-source evidence into a structured context for LLM."""
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from collectors import (
    pubmed, opentargets, intact, uniprot, gwas, chembl, toxicity,
    gnomad, gtex, hpa, dgidb, clinicaltrials, alphafold, reactome,
)

MAX_RETRIES = 3          # 最大リトライ回数
RETRY_WAIT  = [2, 5, 10]  # 待機秒数（指数バックオフ）


def _run_with_retry(fn, key: str, max_retries: int, log):
    """fn を最大 max_retries 回リトライして実行する。"""
    last_err = None
    for attempt in range(max_retries):
        try:
            result = fn()
            if attempt > 0:
                log(f"{key}: OK (リトライ {attempt} 回目で成功)")
            else:
                log(f"{key}: OK")
            return result, None
        except Exception as e:
            last_err = e
            wait = RETRY_WAIT[min(attempt, len(RETRY_WAIT) - 1)]
            log(f"{key}: エラー (試行 {attempt + 1}/{max_retries}) — {e}  → {wait}秒後に再試行")
            if attempt < max_retries - 1:
                time.sleep(wait)
    return None, str(last_err)


def _result_summary(key: str, result, err: str) -> tuple:
    """Returns (status_char, detail_string) for post-collection display."""
    if err:
        short_err = err.split('\n')[0][:55]
        return "✗", short_err
    if result is None:
        return "—", "データなし"
    try:
        if key == "pubmed":
            n = len(result) if isinstance(result, list) else 0
            return "✓", f"{n} 件"
        if key == "opentargets":
            score = result.get("association_score")
            n_drugs = len(result.get("known_drugs", []))
            s = f"score={score:.3f}" if score is not None else "score=N/A"
            return "✓", f"{s}  known_drugs={n_drugs}"
        if key == "uniprot":
            uid  = result.get("uniprot_id", "")
            name = (result.get("protein_name") or "")[:45]
            return "✓", f"{uid}  {name}"
        if key == "intact":
            n = len(result) if isinstance(result, list) else 0
            return "✓", f"{n} interactions"
        if key == "gwas":
            n = len(result) if isinstance(result, list) else 0
            return "✓", f"{n} hits"
        if key == "clinvar":
            n = len(result) if isinstance(result, list) else 0
            return "✓", f"{n} variants"
        if key == "chembl":
            n = len(result) if isinstance(result, list) else 0
            return "✓", f"{n} drugs"
        if key == "gnomad":
            pli   = result.get("pLI")
            loeuf = result.get("LOEUF")
            ess   = result.get("essentiality", "")
            pli_s   = f"{float(pli):.3f}"   if pli   is not None else "N/A"
            loeuf_s = f"{float(loeuf):.3f}" if loeuf is not None else "N/A"
            return "✓", f"pLI={pli_s}  LOEUF={loeuf_s}  {ess}"
        if key == "gtex":
            top = result.get("top_tissues", [])
            s = f"top: {top[0]['tissue']}({top[0]['tpm']:.0f} TPM)" if top else "no data"
            return "✓", s
        if key == "hpa":
            n      = len(result.get("tissue_expression", []))
            subcell = ", ".join(result.get("subcellular", [])[:2]) or "N/A"
            return "✓", f"{n} tissues  loc={subcell}"
        if key == "dgidb":
            items    = result if isinstance(result, list) else []
            approved = sum(1 for d in items if d.get("approved"))
            return "✓", f"{len(items)} interactions ({approved} approved)"
        if key == "clinicaltrials":
            n = len(result) if isinstance(result, list) else 0
            return "✓", f"{n} trials"
        if key == "alphafold":
            plddt = result.get("mean_plddt")
            conf  = result.get("confidence", "")
            return "✓", f"pLDDT={plddt}  {conf}"
        if key == "reactome":
            items   = result if isinstance(result, list) else []
            disease = sum(1 for p in items if p.get("is_disease"))
            return "✓", f"{len(items)} pathways ({disease} disease-related)"
        if key == "toxicity":
            pb      = (result or {}).get("pubchem_bioassay", {})
            ae      = (result or {}).get("drug_adverse_events", {})
            return "✓", f"assays={pb.get('assay_count', 0)}  AE drugs={len(ae)}"
    except Exception:
        pass
    return "✓", ""


def _trunc_at_sentence(text: str, max_chars: int) -> str:
    """Truncate text at the last sentence boundary within max_chars.

    Falls back to word boundary, then hard cut, to avoid mid-word breaks.
    """
    if not text or len(text) <= max_chars:
        return text
    window = text[:max_chars]
    # prefer ending after '. ', '.\n', '! ', '? '
    for sep in ('. ', '.\n', '! ', '? '):
        idx = window.rfind(sep)
        if idx >= int(max_chars * 0.55):
            return window[:idx + 1]
    # fall back to word boundary
    idx = window.rfind(' ')
    if idx >= int(max_chars * 0.55):
        return window[:idx] + ' ...'
    return window + ' ...'


def collect_all(
    gene: str,
    disease: str,
    verbose: bool = True,
    gene_id: str = None,
    disease_id: str = None,
    max_retries: int = MAX_RETRIES,
) -> dict:
    """Run all collectors in parallel with retry, and return aggregated evidence.

    gene_id / disease_id: optional Ensembl / EFO IDs from OpenTargets widget.
    When provided, they skip redundant ID-resolution queries inside opentargets.py.
    max_retries: 各コレクターの最大リトライ回数（デフォルト 3 回）
    """

    def log(msg):
        if verbose:
            print(f"  [+] {msg}")

    tasks = {
        "pubmed":          lambda: pubmed.search_pubmed(gene, disease, max_results=8),
        "opentargets":     lambda: opentargets.get_target_disease_evidence(
                               gene, disease, gene_id=gene_id, disease_id=disease_id),
        "uniprot":         lambda: uniprot.get_protein_info(gene),
        "intact":          lambda: intact.get_interactions(gene, max_results=15),
        "gwas":            lambda: gwas.get_gwas_associations(gene, disease),
        "clinvar":         lambda: gwas.get_clinvar_variants(gene),
        "chembl":          lambda: chembl.get_drugs_for_target(gene),
        "gnomad":          lambda: gnomad.get_constraint(gene),
        "gtex":            lambda: gtex.get_tissue_expression(gene),
        "hpa":             lambda: hpa.get_expression_profile(gene),
        "dgidb":           lambda: dgidb.get_interactions(gene),
        "clinicaltrials":  lambda: clinicaltrials.get_trials(gene, disease),
        "alphafold":       lambda: alphafold.get_structure_info(gene),
        "reactome":        lambda: reactome.get_pathways(gene),
    }

    results = {}
    errors  = {}

    def _task(key, fn):
        return key, *_run_with_retry(fn, key, max_retries, log)

    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(_task, k, fn): k for k, fn in tasks.items()}
        for future in as_completed(futures):
            key, result, err = future.result()
            results[key] = result
            if err:
                errors[key] = err
                log(f"{key}: 最終失敗 — {err}")

    # Toxicity requires known drugs from chembl/opentargets
    known_drugs = []
    if results.get("chembl"):
        known_drugs.extend(results["chembl"])
    if results.get("opentargets") and isinstance(results["opentargets"], dict):
        known_drugs.extend(results["opentargets"].get("known_drugs", []))

    tox_result, tox_err = _run_with_retry(
        lambda: toxicity.assess_target_safety(gene, known_drugs),
        "toxicity", max_retries, log,
    )
    results["toxicity"] = tox_result
    if tox_err:
        errors["toxicity"] = tox_err

    # ── 収集結果サマリー ─────────────────────────────────────────────────────
    if verbose:
        ORDER = [
            "pubmed", "opentargets", "uniprot", "gwas", "clinvar",
            "chembl", "intact", "gnomad", "gtex", "hpa", "dgidb",
            "clinicaltrials", "alphafold", "reactome", "toxicity",
        ]
        w = 55
        print(f"\n  {'─'*w}")
        print(f"  {'ソース':<16} {'状態':<3} 内容")
        print(f"  {'─'*w}")
        for k in ORDER:
            st, detail = _result_summary(k, results.get(k), errors.get(k))
            icon = "✓" if st == "✓" else ("✗" if st == "✗" else "—")
            print(f"  {k:<16} {icon:<3} {detail}")
        n_ok  = sum(1 for k in ORDER if k not in errors and results.get(k) is not None)
        n_err = len(errors)
        print(f"  {'─'*w}")
        print(f"  完了: {n_ok}/{len(ORDER)} ソース取得成功"
              + (f"  ⚠ エラー: {n_err}件" if n_err else ""))

    return {
        "gene": gene,
        "disease": disease,
        "evidence": results,
        "collection_errors": errors,
    }


# デフォルトのコンテキスト設定
DEFAULT_CONTEXT_CONFIG = dict(
    max_papers       = 5,    # 論文数（PubMed）
    abstract_chars   = 600,  # アブストラクト1件あたりの文字数
    max_drugs        = 8,    # 薬剤数（ChEMBL+OpenTargets合計）
    max_gwas         = 5,    # GWAS ヒット数
    max_clinvar      = 5,    # ClinVar バリアント数
    max_interactions = 10,   # PPI インタラクター数
    max_trials       = 6,    # 臨床試験数
    max_reactome     = 10,   # Reactome パスウェイ数
    gtex_top_n       = 5,    # GTEx 上位組織数
    hpa_top_n        = 8,    # HPA 組織数
    max_dgidb        = 8,    # DGIdb 薬剤-遺伝子相互作用数
    uniprot_chars    = 500,  # UniProt function 文字数
    uniprot_keywords = 10,   # UniProt キーワード数
    uniprot_go_terms = 8,    # UniProt GO term 数
)


def build_llm_context(aggregated: dict, config: dict = None) -> str:
    """Convert aggregated evidence into compact structured context for LLM.

    References use short inline format to keep token count low.
    Full citations are stored separately in aggregated["full_references"].

    config: dict to override DEFAULT_CONTEXT_CONFIG values.
    """
    cfg = {**DEFAULT_CONTEXT_CONFIG, **(config or {})}

    gene = aggregated["gene"]
    disease = aggregated["disease"]
    ev = aggregated["evidence"]

    sections = [
        f"# Evidence: {gene} × {disease}\n"
        f"(Cite sources inline using tags like [Paper 1], [ClinVar 2], [UniProt], etc.)\n"
    ]

    # ── Reference registry ─────────────────────────────────────────────────
    # category -> [(tag, short_cite, full_cite), ...]
    ref_reg: dict[str, list] = {"paper": [], "disease": [], "gene": [], "drug": []}
    ref_cnt: dict[str, int]  = {k: 0 for k in ref_reg}

    def add_ref(category: str, prefix: str, short: str, full: str) -> str:
        """Register citation, return tag. short = 1-line; full = Vancouver."""
        if prefix:
            ref_cnt[category] += 1
            tag = f"[{prefix} {ref_cnt[category]}]"
        else:
            tag = f"[{category}]"
        ref_reg[category].append((tag, short, full))
        return tag

    def _short(author_list, year, journal, pmid="", url="") -> str:
        first = author_list[0] if author_list else "Unknown"
        et_al = " et al." if len(author_list) > 1 else ""
        pmid_str = f" PMID:{pmid}" if pmid else ""
        return f"{first}{et_al} {journal} {year}.{pmid_str} {url}".strip()

    # ── Gene/Protein info ──────────────────────────────────────────────────
    uni = ev.get("uniprot") or {}
    if uni and "error" not in uni:
        uid = uni.get("uniprot_id", "")
        url = f"https://www.uniprot.org/uniprotkb/{uid}" if uid else "https://www.uniprot.org/"
        short = f"UniProt entry {gene} ({uid}). {url}"
        full  = (f"UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. "
                 f"Nucleic Acids Res. 2023;51(D1):D523-D531. Entry: {gene} ({uid}). {url}")
        ref = add_ref("gene", "UniProt", short, full)
        func = _trunc_at_sentence(uni.get("function") or "", cfg["uniprot_chars"])
        kws  = ", ".join(uni.get("keywords", [])[:cfg["uniprot_keywords"]]) or "N/A"
        gos  = "; ".join(g["term"] for g in uni.get("go_terms", [])[:cfg["uniprot_go_terms"]]) or "N/A"
        sections.append(
            f"## Gene/Protein {ref}\n"
            f"- Name: {uni.get('protein_name', 'N/A')}\n"
            f"- Function: {func}\n"
            f"- Location: {', '.join(uni.get('subcellular_location', [])[:4]) or 'N/A'}\n"
            f"- Keywords: {kws}\n"
            f"- GO terms: {gos}\n"
        )

    # ── OpenTargets ────────────────────────────────────────────────────────
    ot = ev.get("opentargets") or {}
    if ot and "error" not in ot:
        ensg = ot.get("ensembl_id", "")
        efo  = ot.get("disease_id", "")
        url  = (f"https://platform.opentargets.org/evidence/{ensg}/{efo}"
                if ensg and efo else "https://platform.opentargets.org/")
        short = f"OpenTargets {gene}×{disease}. {url}"
        full  = (f"Ochoa D, et al. Open Targets Platform. "
                 f"Nucleic Acids Res. 2021;49(D1):D1302-D1310. {url}")
        ref = add_ref("disease", "OpenTargets", short, full)
        score = ot.get("association_score")
        score_str = f"{score:.3f}" if score is not None else "N/A"
        dt_str = " | ".join(f"{k}:{v:.2f}" for k, v in (ot.get("datatype_scores") or {}).items())
        sections.append(
            f"## OpenTargets {ref}\n"
            f"- Score: {score_str} | {dt_str}\n"
        )

    # ── GWAS ───────────────────────────────────────────────────────────────
    gwas_hits = ev.get("gwas") or []
    if gwas_hits:
        lines = []
        for h in gwas_hits[:cfg["max_gwas"]]:
            sid  = h.get("study_id", "")
            year = (h.get("pub_date") or "")[:4]
            auth = h.get("first_author", "")
            url  = f"https://www.ebi.ac.uk/gwas/studies/{sid}" if sid else "https://www.ebi.ac.uk/gwas/"
            short = f"{auth} et al. GWAS Catalog {sid} {year}. {url}"
            full  = (f"{auth} et al. GWAS Catalog Study {sid}: {h.get('trait','')}. {year}. "
                     f"Buniello A et al. Nucleic Acids Res. 2019;47(D1):D1005-D1012. {url}")
            ref = add_ref("disease", "GWAS", short, full)
            snps = ", ".join(h.get("snps", [])[:2])
            lines.append(f"  - {h['trait']} p={h['p_value']} OR={h['or_beta']} SNPs:{snps} {ref}")
        sections.append(f"## GWAS\n" + "\n".join(lines) + "\n")

    # ── ClinVar ────────────────────────────────────────────────────────────
    cv_hits = ev.get("clinvar") or []
    if cv_hits:
        lines = []
        for v in cv_hits[:cfg["max_clinvar"]]:
            vid = str(v.get("uid", ""))
            url = (f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/"
                   if vid else f"https://www.ncbi.nlm.nih.gov/clinvar/?term={gene}")
            short = f"ClinVar VarID:{vid} {v.get('clinical_significance','')}. {url}"
            full  = (f"Landrum MJ et al. ClinVar. Nucleic Acids Res. 2020;48(D1):D835-D844. "
                     f"Variant: {v.get('title','')[:80]} [VarID:{vid}]. {url}")
            ref = add_ref("disease", "ClinVar", short, full)
            lines.append(f"  - {v['title'][:65]} | {v['clinical_significance']} | {v['condition']} {ref}")
        sections.append(f"## ClinVar Variants\n" + "\n".join(lines) + "\n")

    # ── Existing drugs ─────────────────────────────────────────────────────
    drugs_chembl = ev.get("chembl") or []
    drugs_ot     = ot.get("known_drugs", []) if ot else []
    all_drugs = {
        d.get("name") or d.get("drug", ""): d
        for d in drugs_chembl + drugs_ot
        if d.get("name") or d.get("drug")
    }
    if all_drugs:
        lines = []
        for name, d in list(all_drugs.items())[:cfg["max_drugs"]]:
            cid   = d.get("chembl_id", "")
            url   = (f"https://www.ebi.ac.uk/chembl/compound_report_card/{cid}/"
                     if cid else "https://www.ebi.ac.uk/chembl/")
            phase = d.get("max_phase") or d.get("phase")
            mech  = d.get("mechanism") or d.get("mechanism_of_action", "")
            short = f"{name} phase:{phase} {mech[:60]}. {url}"
            full  = (f"Mendez D et al. ChEMBL. Nucleic Acids Res. 2019;47(D1):D930-D940. "
                     f"Compound: {name}" + (f" ({cid})" if cid else "")
                     + (f" phase:{phase}" if phase else "")
                     + (f" mech:{mech[:80]}" if mech else "") + f". {url}")
            ref = add_ref("drug", "ChEMBL", short, full)
            lines.append(f"  - {name} | Ph:{phase} | {mech[:60] or 'N/A'} {ref}")
        sections.append(f"## Drugs targeting {gene}\n" + "\n".join(lines) + "\n")
    else:
        sections.append(f"## Drugs\nNo approved drugs found for {gene}.\n")

    # ── Protein interactions ───────────────────────────────────────────────
    interactions = ev.get("intact") or []
    if interactions:
        url   = f"https://www.ebi.ac.uk/intact/search?query={gene}"
        short = f"IntAct PPI for {gene}. {url}"
        full  = (f"Orchard S et al. IntAct. Nucleic Acids Res. 2014;42(D1):D358-D363. {url}")
        ref   = add_ref("gene", "IntAct", short, full)
        partners = list(dict.fromkeys(
            p for ix in interactions[:cfg["max_interactions"]] for p in ix.get("partners", [])
        ))[:cfg["max_interactions"]]
        sections.append(
            f"## PPI (IntAct) {ref}\n"
            f"- Interactors: {', '.join(partners) or 'N/A'}\n"
        )

    # ── Toxicity ───────────────────────────────────────────────────────────
    tox = ev.get("toxicity") or {}
    pb  = tox.get("pubchem_bioassay", {})
    ae  = tox.get("drug_adverse_events", {})
    url = f"https://pubchem.ncbi.nlm.nih.gov/#query={gene}&input_type=gene"
    short = f"PubChem BioAssay {gene}. {url}"
    full  = f"Kim S et al. PubChem 2023. Nucleic Acids Res. 2023;51(D1):D1373-D1380. {url}"
    ref_pc = add_ref("gene", "PubChem", short, full)
    ae_str = "; ".join(
        f"{drug}: " + ", ".join(f"{e['reaction']}({e['count']})" for e in evts[:2])
        for drug, evts in list(ae.items())[:2]
    ) or "N/A"
    sections.append(
        f"## Safety {ref_pc}\n"
        f"- BioAssays: {pb.get('assay_count', 0)} | AEs: {ae_str}\n"
    )

    # ── Literature ─────────────────────────────────────────────────────────
    papers = ev.get("pubmed") or []
    if papers:
        paper_blocks = []
        for p in papers[:cfg["max_papers"]]:
            pmid     = p.get("pmid", "")
            title    = p.get("title", "")
            journal  = p.get("journal", "")
            year     = p.get("year", "")
            authors  = p.get("authors", [])
            abstract = (p.get("abstract") or "").strip()
            url      = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
            first3   = authors[:3]
            auth_str = ", ".join(first3) + (", et al" if len(authors) > 3 else "")
            short    = _short(authors, year, journal, pmid, url)
            full     = f"{auth_str}. {title} {journal}. {year}. PMID:{pmid}. {url}"
            ref = add_ref("paper", "Paper", short, full)
            snippet = _trunc_at_sentence(abstract, cfg["abstract_chars"]) \
                      if abstract else "(no abstract)"
            paper_blocks.append(
                f"### {ref} {title[:80]} ({year})\n"
                f"_{auth_str} | {journal}_\n"
                f"{snippet}\n"
            )
        sections.append("## Literature\n\n" + "\n".join(paper_blocks))

    # ── gnomAD ────────────────────────────────────────────────────────────
    gnom = ev.get("gnomad") or {}
    if gnom and "error" not in gnom:
        url   = gnom.get("url", "https://gnomad.broadinstitute.org/")
        short = f"gnomAD {gene} pLI={gnom.get('pLI')} LOEUF={gnom.get('LOEUF')}. {url}"
        full  = (f"Chen S et al. Nature. 2024;625:92-100. {short}")
        ref   = add_ref("gene", "gnomAD", short, full)
        sections.append(
            f"## Constraint (gnomAD) {ref}\n"
            f"- pLI:{gnom.get('pLI')} LOEUF:{gnom.get('LOEUF')} "
            f"missenseOE:{gnom.get('oe_missense')} — {gnom.get('essentiality')}\n"
        )

    # ── GTEx ──────────────────────────────────────────────────────────────
    gtex_data = ev.get("gtex") or {}
    if gtex_data and "error" not in gtex_data:
        url   = gtex_data.get("url", "https://gtexportal.org/")
        short = f"GTEx expression {gene}. {url}"
        full  = f"GTEx Consortium. Science. 2020;369(6509):1318-1330. {short}"
        ref   = add_ref("gene", "GTEx", short, full)
        top3  = ", ".join(
            f"{t['tissue']}({t['tpm']:.0f})" for t in gtex_data.get("top_tissues", [])[:cfg["gtex_top_n"]]
        )
        key   = "; ".join(
            f"{t['tissue']}:{t['tpm']:.0f}" for t in gtex_data.get("key_tissues", []) if t["tpm"] > 0
        ) or "N/A"
        sections.append(
            f"## Expression GTEx {ref}\n"
            f"- Top: {top3} | Safety tissues: {key}\n"
        )

    # ── Human Protein Atlas ────────────────────────────────────────────────
    hpa_data = ev.get("hpa") or {}
    if hpa_data and "error" not in hpa_data:
        url   = hpa_data.get("url", "https://www.proteinatlas.org/")
        short = f"HPA {gene}. {url}"
        full  = f"Uhlén M et al. Science. 2015;347(6220):1260419. {short}"
        ref   = add_ref("gene", "HPA", short, full)
        subcell    = ", ".join(hpa_data.get("subcellular", [])[:4]) or "N/A"
        prot_class = ", ".join(hpa_data.get("protein_class", [])[:3]) or "N/A"
        tissues    = " | ".join(
            f"{t['tissue']}:{t['level']}"
            for t in hpa_data.get("tissue_expression", [])[:cfg["hpa_top_n"]]
        ) or "N/A"
        sections.append(
            f"## Protein Atlas {ref}\n"
            f"- Class:{prot_class} | Location:{subcell}\n"
            f"- Expression: {tissues}\n"
        )

    # ── DGIdb ─────────────────────────────────────────────────────────────
    dgi_data = ev.get("dgidb") or []
    if dgi_data:
        url   = f"https://dgidb.org/genes/{gene}#interactions"
        short = f"DGIdb {gene}. {url}"
        full  = f"Cannon M et al. Nucleic Acids Res. 2024;52(D1):D1227-D1235. {short}"
        ref   = add_ref("drug", "DGIdb", short, full)
        approved = [d for d in dgi_data if d.get("approved")]
        rows = " | ".join(
            f"{d['drug_name']}({'✓' if d.get('approved') else 'inv'},{d.get('interaction_type','')})"
            for d in dgi_data[:cfg["max_dgidb"]]
        )
        sections.append(
            f"## DGIdb {ref}\n"
            f"- {len(dgi_data)} interactions ({len(approved)} approved): {rows}\n"
        )

    # ── ClinicalTrials ────────────────────────────────────────────────────
    ct_data = ev.get("clinicaltrials") or []
    if ct_data:
        url   = f"https://clinicaltrials.gov/search?cond={disease.replace(' ','%20')}&intr={gene}"
        short = f"ClinicalTrials {gene}×{disease}. {url}"
        full  = f"ClinicalTrials.gov. U.S. NLM. {short}"
        ref   = add_ref("disease", "ClinicalTrials", short, full)
        rows  = " | ".join(
            f"{t['nct_id']}({t['phase']},{t['status'][:8]})"
            for t in ct_data[:cfg["max_trials"]]
        )
        sections.append(
            f"## Clinical Trials {ref}\n"
            f"- {len(ct_data)} trials: {rows}\n"
        )

    # ── AlphaFold ─────────────────────────────────────────────────────────
    af_data = ev.get("alphafold") or {}
    if af_data and "error" not in af_data:
        url   = af_data.get("view_url", "https://alphafold.ebi.ac.uk/")
        short = f"AlphaFold {af_data.get('entry_id','')} pLDDT={af_data.get('mean_plddt')}. {url}"
        full  = f"Jumper J et al. Nature. 2021;596:583-589. {short}"
        ref   = add_ref("gene", "AlphaFold", short, full)
        sections.append(
            f"## Structure AlphaFold {ref}\n"
            f"- pLDDT:{af_data.get('mean_plddt')} — {af_data.get('confidence')}\n"
        )

    # ── Reactome ──────────────────────────────────────────────────────────
    react_data = ev.get("reactome") or []
    if react_data:
        url   = f"https://reactome.org/content/query?q={gene}&species=Homo+sapiens"
        short = f"Reactome pathways {gene}. {url}"
        full  = f"Milacic M et al. Nucleic Acids Res. 2024;52(D1):D672-D678. {short}"
        ref   = add_ref("gene", "Reactome", short, full)
        dpw   = [p for p in react_data if p.get("is_disease")]
        names = " | ".join(p["name"] for p in react_data[:cfg["max_reactome"]])
        sections.append(
            f"## Reactome {ref}\n"
            f"- {len(react_data)} pathways ({len(dpw)} disease): {names}\n"
        )

    # ── References (compact) ──────────────────────────────────────────────
    HEADERS = {
        "paper":   "### Papers",
        "disease": "### Disease DBs (GWAS/ClinVar/OpenTargets/ClinicalTrials)",
        "gene":    "### Gene/Protein DBs (UniProt/IntAct/gnomAD/GTEx/HPA/AlphaFold/Reactome/PubChem)",
        "drug":    "### Drug DBs (ChEMBL/DGIdb)",
    }
    ref_lines = ["## References\n"]
    for cat, header in HEADERS.items():
        entries = ref_reg.get(cat, [])
        if entries:
            ref_lines.append(header)
            for tag, short, _ in entries:
                ref_lines.append(f"{tag} {short}")
            ref_lines.append("")

    sections.append("\n".join(ref_lines))

    # 完全引用を別フィールドに保存（レポート保存時に使用）
    aggregated["full_references"] = {
        cat: [(tag, full) for tag, _, full in entries]
        for cat, entries in ref_reg.items()
    }

    return "\n".join(sections)
