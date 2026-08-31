"""Aggregate multi-source evidence into a structured context for LLM."""
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from collectors import (
    pubtator, opentargets, intact, uniprot, gwas, chembl, toxicity,
    gnomad, gtex, hpa, dgidb, clinicaltrials, alphafold, reactome, gprofiler,
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
            n_clin = sum(1 for p in result if p.get("is_clinical")) if isinstance(result, list) else 0
            return "✓", f"{n} 件 (臨床 {n_clin} 件)"
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
            if "error" in result:
                return "✗", result["error"][:55]
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
            tc      = (result or {}).get("toxcast", {})
            ae      = (result or {}).get("drug_adverse_events", {})
            return "✓", f"toxcast_assays={tc.get('assay_count', 0)}  AE drugs={len(ae)}"
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
        "pubmed":          lambda: pubtator.search_pubtator(
                               gene, disease, max_results=100,
                               disease_efo_id=disease_id),
        "opentargets":     lambda: opentargets.get_target_disease_evidence(
                               gene, disease, gene_id=gene_id, disease_id=disease_id),
        "uniprot":         lambda: uniprot.get_protein_info(gene),
        "intact":          lambda: intact.get_interactions(gene, max_results=15),
        "gwas":            lambda: gwas.get_gwas_associations(gene, disease),
        "chembl":          lambda: chembl.get_drugs_for_target(gene),
        "gnomad":          lambda: gnomad.get_constraint(gene),
        "gtex":            lambda: gtex.get_tissue_expression(gene),
        "hpa":             lambda: hpa.get_expression_profile(gene),
        "dgidb":           lambda: dgidb.get_interactions(gene),
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

    # ClinVar: opentargets の synonyms を使い疾患関連バリアントに絞り込む
    ot_result = results.get("opentargets")
    disease_synonyms = (
        ot_result.get("disease_synonyms", [])
        if isinstance(ot_result, dict) else []
    )
    cv_result, cv_err = _run_with_retry(
        lambda: gwas.get_clinvar_variants(gene, disease_query=disease, disease_synonyms=disease_synonyms),
        "clinvar", max_retries, log,
    )
    results["clinvar"] = cv_result
    if cv_err:
        errors["clinvar"] = cv_err

    # ClinicalTrials: gene-symbol 検索に加え、known_drugs の薬剤名でも検索する
    # （治験の大半は標的遺伝子名ではなく薬剤コード名でしか言及されないため、
    #  遺伝子名単独の検索では大半の治験を取りこぼす）
    drug_names = [d.get("name") or d.get("drug") or "" for d in known_drugs]
    ct_result, ct_err = _run_with_retry(
        lambda: clinicaltrials.get_trials(gene, disease, drug_names=drug_names,
                                           disease_efo_id=disease_id),
        "clinicaltrials", max_retries, log,
    )
    results["clinicaltrials"] = ct_result
    if ct_err:
        errors["clinicaltrials"] = ct_err

    # パスウェイ経由の関連解析（opentargets 上位遺伝子 × Reactome 共有パスウェイ）
    ot_disease_id = ot_result.get("disease_id") if isinstance(ot_result, dict) else None
    if ot_disease_id:
        try:
            disease_genes = opentargets.get_disease_top_genes(ot_disease_id, top_n=20)
            pathway_connections = reactome.find_pathway_connections(
                gene, disease_genes, max_partners=5
            )
            results["pathway_connections"] = pathway_connections
            log(f"pathway_connections: {len(pathway_connections)} 件")
        except Exception as e:
            errors["pathway_connections"] = str(e)
            results["pathway_connections"] = []
            log(f"pathway_connections: FAILED ({e})")
    else:
        results["pathway_connections"] = []

    # Disease pathway enrichment + target gene pathway fit assessment
    if ot_disease_id:
        try:
            disease_genes_for_enrich = opentargets.get_disease_top_genes(ot_disease_id, top_n=20)
            enriched = gprofiler.enrich_gene_list([g["symbol"] for g in disease_genes_for_enrich])
            disease_pathway_ids = {p["term_id"] for p in enriched if p["source"] == "REAC"}
            target_in = reactome.get_gene_pathway_membership(gene, disease_pathway_ids)
            score = len(target_in) / max(1, min(20, len(disease_pathway_ids)))
            results["pathway_fit"] = {
                "disease_pathways": enriched[:20],
                "target_in_disease_pathways": target_in,
                "pathway_overlap_score": round(score, 3),
                "gene_list_size": len(disease_genes_for_enrich),
            }
            log(f"pathway_fit: score={results['pathway_fit']['pathway_overlap_score']}, "
                f"target_in={len(target_in)} pathways")
        except Exception as e:
            errors["pathway_fit"] = str(e)
            results["pathway_fit"] = {
                "disease_pathways": [], "target_in_disease_pathways": [],
                "pathway_overlap_score": 0.0, "gene_list_size": 0,
            }
            log(f"pathway_fit: FAILED ({e})")
    else:
        results["pathway_fit"] = {
            "disease_pathways": [], "target_in_disease_pathways": [],
            "pathway_overlap_score": 0.0, "gene_list_size": 0,
        }

    # パスウェイ隣接遺伝子の論文を補足取得（並列実行）
    partners = [
        conn.get("partner", "")
        for conn in (results.get("pathway_connections") or [])[:3]
        if conn.get("partner")
    ]

    def _fetch_partner_papers(partner: str):
        try:
            papers = pubtator.search_pubtator(
                partner, disease, max_results=5, disease_efo_id=disease_id
            )
            papers = [p for p in papers if p.get("relevance_score", 0) > 0]
            return partner, papers
        except Exception as e:
            log(f"related papers ({partner}): FAILED ({e})")
            return partner, []

    partner_papers: dict[str, list] = {}
    if partners:
        with ThreadPoolExecutor(max_workers=len(partners)) as ex:
            for partner, papers in ex.map(_fetch_partner_papers, partners):
                if papers:
                    partner_papers[partner] = papers
                    log(f"related papers ({partner}): {len(papers)} 件")
    results["related_gene_papers"] = partner_papers

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
    max_papers       = 8,    # 論文数（PubMed、臨床妥当性/MoA妥当性に均等配分）
    abstract_chars   = 600,  # アブストラクト1件あたりの文字数
    max_drugs        = 8,    # 薬剤数（ChEMBL+OpenTargets合計）
    max_gwas         = 5,    # GWAS ヒット数
    max_clinvar      = 20,   # ClinVar バリアント数
    max_interactions = 10,   # PPI インタラクター数
    max_trials       = 10,   # 臨床試験数（context に含める直近件数。全件数自体は別途表示）
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

    def add_ref(category: str, prefix: str, short: str, full: str, url: str = "") -> str:
        """Register citation, return tag. short = 1-line; full = Vancouver.

        url: リンク可能な一次情報源のURL。References セクションはこれを使い
             サーバー側で確定的に生成する（LLM 生成には任せない）。
        """
        if prefix:
            ref_cnt[category] += 1
            tag = f"[{prefix} {ref_cnt[category]}]"
        else:
            tag = f"[{category}]"
        ref_reg[category].append((tag, short, full, url))
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
        ref = add_ref("gene", "UniProt", short, full, url)
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
        ref = add_ref("disease", "OpenTargets", short, full, url)
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
            # collectors/gwas.py が返す個々の association には study_id/pub_date/
            # first_author は含まれない（GWAS Catalog API がSNPアソシエーション単位
            # では返さないため）。実際にリンク先として使えるのは efoTraits 経由で
            # 取得した trait ページ URL (gwas_url) のみ。
            url  = h.get("gwas_url") or "https://www.ebi.ac.uk/gwas/"
            short = f"GWAS Catalog: {h.get('trait','')}. {url}"
            full  = (f"GWAS Catalog — {h.get('trait','')}. "
                     f"Buniello A et al. Nucleic Acids Res. 2019;47(D1):D1005-D1012. {url}")
            ref = add_ref("disease", "GWAS", short, full, url)
            snps = ", ".join(h.get("snps", [])[:2])
            lines.append(f"  - {h['trait']} p={h['p_value']} OR={h['or_beta']} SNPs:{snps} {ref}")
        sections.append(f"## GWAS\n" + "\n".join(lines) + "\n")

    # ── ClinVar ────────────────────────────────────────────────────────────
    cv_hits = ev.get("clinvar") or []
    if cv_hits:
        lines = []
        for v in cv_hits[:cfg["max_clinvar"]]:
            vid = str(v.get("variant_id") or v.get("uid", ""))
            url = (f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/"
                   if vid else f"https://www.ncbi.nlm.nih.gov/clinvar/?term={gene}")
            short = f"ClinVar VarID:{vid} {v.get('clinical_significance','')}. {url}"
            full  = (f"Landrum MJ et al. ClinVar. Nucleic Acids Res. 2020;48(D1):D835-D844. "
                     f"Variant: {v.get('title','')[:80]} [VarID:{vid}]. {url}")
            ref = add_ref("disease", "ClinVar", short, full, url)
            lines.append(
                f"  - {v.get('title','')[:65]} | {v.get('clinical_significance','')} "
                f"| {v.get('condition','')} | {v.get('review_status','')} {ref}"
            )
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
            ref = add_ref("drug", "ChEMBL", short, full, url)
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
        ref   = add_ref("gene", "IntAct", short, full, url)
        partners = list(dict.fromkeys(
            p for ix in interactions[:cfg["max_interactions"]] for p in ix.get("partners", [])
        ))[:cfg["max_interactions"]]
        sections.append(
            f"## PPI (IntAct) {ref}\n"
            f"- Interactors: {', '.join(partners) or 'N/A'}\n"
        )

    # ── Toxicity ───────────────────────────────────────────────────────────
    tox = ev.get("toxicity") or {}
    tc  = tox.get("toxcast", {})
    ae  = tox.get("drug_adverse_events", {})
    url = "https://comptox.epa.gov/dashboard"
    short = f"EPA ToxCast {gene}. {url}"
    full  = f"US EPA. ToxCast/Tox21 Bioactivity Data via CTX Bioactivity API. {url}"
    ref_tc = add_ref("gene", "ToxCast", short, full, url)
    ae_str = "; ".join(
        f"{drug}: " + ", ".join(f"{e['reaction']}({e['count']})" for e in evts[:2])
        for drug, evts in list(ae.items())[:2]
    ) or "N/A"
    tc_str = f"{tc.get('assay_count', 0)}" if tc.get("available") else "N/A (no API key)"
    sections.append(
        f"## Safety {ref_tc}\n"
        f"- ToxCast Assays: {tc_str} | AEs: {ae_str}\n"
    )

    # ── Literature ─────────────────────────────────────────────────────────
    # 論文は2つの妥当性軸で使う:
    #   1. 臨床妥当性 (Clinical Validity)  — 実際にヒトで臨床試験が行われた実例
    #      があるか（is_clinical=True の論文）
    #   2. MoA妥当性 (Mechanistic Validity) — 標的と疾患を結びつける機序・機能
    #      的根拠があるか（前臨床/機序論文, is_clinical=False）
    # どちらか一方だけに偏らないよう、max_papers の枠を両カテゴリに配分する
    # （臨床論文を優先しすぎると機序の説明が抜け落ちるため）。
    papers = ev.get("pubmed") or []
    if papers:
        def _paper_block(p):
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
            ref = add_ref("paper", "Paper", short, full, url)
            snippet = _trunc_at_sentence(abstract, cfg["abstract_chars"]) \
                      if abstract else "(no abstract)"
            return f"### {ref} {title[:80]} ({year})\n_{auth_str} | {journal}_\n{snippet}\n"

        clinical_papers = [p for p in papers if p.get("is_clinical")]
        moa_papers      = [p for p in papers if not p.get("is_clinical")]
        half = max(1, cfg["max_papers"] // 2)
        clin_pick = clinical_papers[:half]
        moa_pick  = moa_papers[:cfg["max_papers"] - len(clin_pick)]
        # 片方が枠に満たない場合はもう片方の枠を広げて埋める
        if len(clin_pick) < half:
            moa_pick = moa_papers[:cfg["max_papers"] - len(clin_pick)]
        elif len(moa_pick) < cfg["max_papers"] - half:
            clin_pick = clinical_papers[:cfg["max_papers"] - len(moa_pick)]

        blocks = ["## Literature\n"]
        blocks.append(
            f"_Pool: top {len(papers)} most-recent PubTator3 hits with gene+disease entity match "
            f"({len(clinical_papers)} clinical, {len(moa_papers)} preclinical/mechanistic)._\n"
        )
        blocks.append(
            f"### Clinical Validity — evidence of actual clinical trial/patient precedent "
            f"({len(clin_pick)} shown)\n"
        )
        blocks.append("\n".join(_paper_block(p) for p in clin_pick) if clin_pick
                       else "_No clinical-stage literature found for this gene×disease pair._\n")
        blocks.append(
            f"### MoA Validity — mechanistic/functional rationale linking target to disease "
            f"({len(moa_pick)} shown)\n"
        )
        blocks.append("\n".join(_paper_block(p) for p in moa_pick) if moa_pick
                       else "_No mechanistic literature found beyond the clinical papers above._\n")
        sections.append("\n".join(blocks))

    # ── gnomAD ────────────────────────────────────────────────────────────
    gnom = ev.get("gnomad") or {}
    if gnom and "error" not in gnom:
        url   = gnom.get("url", "https://gnomad.broadinstitute.org/")
        short = f"gnomAD {gene} pLI={gnom.get('pLI')} LOEUF={gnom.get('LOEUF')}. {url}"
        full  = (f"Chen S et al. Nature. 2024;625:92-100. {short}")
        ref   = add_ref("gene", "gnomAD", short, full, url)
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
        ref   = add_ref("gene", "GTEx", short, full, url)
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
        ref   = add_ref("gene", "HPA", short, full, url)
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
        ref   = add_ref("drug", "DGIdb", short, full, url)
        approved = [d for d in dgi_data if d.get("approved")]
        rows = " | ".join(
            f"{d['drug_name']}({'✓' if d.get('approved') else 'inv'},{d.get('interaction_type','')})"
            for d in dgi_data[:cfg["max_dgidb"]]
        )
        sections.append(
            f"## DGIdb {ref}\n"
            f"- {len(dgi_data)} interactions ({len(approved)} approved): {rows}\n"
        )

    # ── ClinicalTrials（競合状況＝新規参入リスクの指標として使用） ──────────────
    # 件数に上限は設けず全件取得しているため、ここでの trial 数は競合の激しさを
    # 直接反映する。LLM には「試験数が多いほど新規参入リスクが高い」という
    # 読み方を明示し、context には直近の試験を優先して渡す（ct_data は
    # collectors/clinicaltrials.py で開始日降順に既にソート済み）。
    ct_data = ev.get("clinicaltrials") or []
    if ct_data:
        url   = f"https://clinicaltrials.gov/search?cond={disease.replace(' ','%20')}&term={gene}"
        short = f"ClinicalTrials {gene}×{disease}. {url}"
        full  = f"ClinicalTrials.gov. U.S. NLM. {short}"
        ref   = add_ref("disease", "ClinicalTrials", short, full, url)
        rows  = " | ".join(
            f"{t['nct_id']}({t['phase']},{t['status'][:8]},{t.get('start_date','') or 'N/A'})"
            for t in ct_data[:cfg["max_trials"]]
        )
        n = len(ct_data)
        risk = "HIGH" if n >= 10 else ("MODERATE" if n >= 3 else "LOW")
        sections.append(
            f"## Competitive Landscape {ref}\n"
            f"- Total trials for this target×disease: {n}  →  new-entrant risk: {risk} "
            f"(more competing trials = harder to differentiate / higher risk)\n"
            f"- Most recent {min(cfg['max_trials'], n)} trials: {rows}\n"
        )

    # ── AlphaFold ─────────────────────────────────────────────────────────
    af_data = ev.get("alphafold") or {}
    if af_data and "error" not in af_data:
        url   = af_data.get("view_url", "https://alphafold.ebi.ac.uk/")
        short = f"AlphaFold {af_data.get('entry_id','')} pLDDT={af_data.get('mean_plddt')}. {url}"
        full  = f"Jumper J et al. Nature. 2021;596:583-589. {short}"
        ref   = add_ref("gene", "AlphaFold", short, full, url)
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
        ref   = add_ref("gene", "Reactome", short, full, url)
        dpw   = [p for p in react_data if p.get("is_disease")]
        names = " | ".join(p["name"] for p in react_data[:cfg["max_reactome"]])
        sections.append(
            f"## Reactome {ref}\n"
            f"- {len(react_data)} pathways ({len(dpw)} disease): {names}\n"
        )

    # ── Pathway-level connections (indirect associations) ────────────────────
    pw_connections = ev.get("pathway_connections") or []
    if pw_connections:
        lines = []
        for conn in pw_connections:
            partner = conn.get("partner", "")
            p_score = conn.get("partner_score")
            p_score_str = f"{p_score:.3f}" if p_score is not None else "N/A"
            shared = conn.get("shared_pathways", [])
            pw_names = " | ".join(p["name"] for p in shared[:3])
            lines.append(
                f"  - {gene} ↔ {partner} (disease assoc score: {p_score_str}): {pw_names}"
            )
        sections.append(
            f"## Pathway-level Indirect Connections\n"
            f"Genes strongly associated with {disease} that share pathways with {gene}:\n"
            + "\n".join(lines) + "\n"
            + f"Note: These pathway connections suggest {gene} may influence {disease} "
              f"indirectly through shared biological mechanisms.\n"
        )

    # ── Disease Pathway Analysis (g:Profiler enrichment + target fit) ─────────
    pathway_fit = ev.get("pathway_fit") or {}
    if pathway_fit:
        disease_pathways = pathway_fit.get("disease_pathways") or []
        target_in = pathway_fit.get("target_in_disease_pathways") or []
        overlap_score = pathway_fit.get("pathway_overlap_score", 0.0)
        gene_list_size = pathway_fit.get("gene_list_size", 0)

        pw_lines = [
            f"## Disease Pathway Analysis\n"
            f"Disease-associated pathways (from g:Profiler enrichment of top {gene_list_size} {disease} genes):"
        ]
        for pw in disease_pathways[:10]:
            src = pw.get("source", "")
            name = pw.get("name", "")
            pval = pw.get("p_value", 1.0)
            isect = pw.get("intersection_size", 0)
            tsize = pw.get("term_size", 0)
            pw_lines.append(f"- [{src}] {name} (p={pval:.3g}, {isect} of {tsize} disease genes)")

        total_dp = len(disease_pathways)
        pw_lines.append(f"\nTarget gene {gene} membership in disease pathways:")
        pw_lines.append(
            f"- Direct member of {len(target_in)}/{total_dp} enriched disease pathways "
            f"(overlap score: {overlap_score:.2f})"
        )
        if target_in:
            names = ", ".join(
                f"{p.get('name', p.get('pathway_id', ''))} [{p.get('source','')}]"
                for p in target_in[:5]
            )
            pw_lines.append(f"- Present in: {names}")
            # 共存する疾患関連遺伝子を列挙（LLMが関係性を説明するための材料）
            co_genes: list[str] = []
            for pw in disease_pathways:
                if any(p.get("name") == pw.get("name") for p in target_in):
                    co_genes.extend(pw.get("genes") or [])
            co_genes_uniq = list(dict.fromkeys(g for g in co_genes if g.upper() != gene.upper()))[:10]
            if co_genes_uniq:
                pw_lines.append(f"- Co-pathway disease genes: {', '.join(co_genes_uniq)}")
        else:
            pw_lines.append(
                f"- Not directly present in top disease pathways "
                f"(indirect connection via shared pathway partners)"
            )
        sections.append("\n".join(pw_lines) + "\n")

    # ── Network × Disease Gene Overlap ──────────────────────────────────────
    net_overlap = ev.get("network_disease_overlap") or {}
    if net_overlap and net_overlap.get("overlap_count", 0) > 0:
        ws  = net_overlap.get("weighted_score", 0.0)
        sr  = net_overlap.get("simple_ratio", 0.0)
        ov  = net_overlap.get("overlap_count", 0)
        dg  = net_overlap.get("disease_gene_count", 0)
        pp  = net_overlap.get("ppi_partner_count", 0)
        top = net_overlap.get("overlapping_genes", [])[:8]
        top_str = ", ".join(
            f"{g['symbol']} (OT={g['score']:.2f})" for g in top if g.get("symbol")
        )
        sections.append(
            f"## Network–Disease Gene Overlap\n"
            f"PPI partners of {gene} that are also top {disease} genes in OpenTargets:\n"
            f"- Weighted overlap score: {ws:.3f} (sum of OT scores of overlapping genes / total)\n"
            f"- Simple overlap: {ov}/{dg} disease genes ({sr*100:.1f}%) found in {gene} PPI network ({pp} partners)\n"
            f"- Top overlapping genes: {top_str}\n"
            f"Interpretation: a higher weighted score indicates that {gene} is more central "
            f"to the {disease} disease gene network via direct protein interactions.\n"
        )

    # ── Related gene papers (pathway-connected partners) ────────────────────
    related_gene_papers = ev.get("related_gene_papers") or {}
    for partner_gene, papers in related_gene_papers.items():
        if not papers:
            continue
        lines = [f"## Related Literature: {partner_gene} × {disease}",
                 f"(Pathway-connected gene; papers included for mechanistic context)\n"]
        for p in papers:
            pmid  = p.get("pmid", "")
            title = p.get("title", "")
            year  = p.get("year", "")
            abstract = p.get("abstract", "") or ""
            url   = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
            short = f"{partner_gene} paper — {title[:80]}... ({year})"
            full  = f"{title} ({year})"
            ref = add_ref("paper", "Paper", short, full, url)
            lines.append(f"{ref} [{year}] {title}")
            if abstract:
                lines.append(f"Abstract: {abstract[:300]}...")
            lines.append("")
        sections.append("\n".join(lines))

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
            for tag, short, _full, _url in entries:
                ref_lines.append(f"{tag} {short}")
            ref_lines.append("")

    sections.append("\n".join(ref_lines))

    # 完全引用+URLを別フィールドに保存。最終的な "## References" セクションは
    # LLM に生成させず、report.references_md() でサーバー側から確定的に
    # 生成しリンクを付与する（LLMの自由記述だとタグの欠落・書式崩れが起きるため）。
    aggregated["full_references"] = {
        cat: [(tag, full, url) for tag, _short, full, url in entries]
        for cat, entries in ref_reg.items()
    }

    return "\n".join(sections)
