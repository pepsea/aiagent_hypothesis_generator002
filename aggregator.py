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

    return {
        "gene": gene,
        "disease": disease,
        "evidence": results,
        "collection_errors": errors,
    }


def build_llm_context(aggregated: dict) -> str:
    """Convert aggregated evidence into a structured text context for LLM.

    Each evidence item is assigned a numbered reference [Ref N] so the LLM
    can cite specific sources in the hypothesis report.
    """
    gene = aggregated["gene"]
    disease = aggregated["disease"]
    ev = aggregated["evidence"]

    sections = [
        f"# Evidence Summary for Drug Target Hypothesis\n"
        f"Gene: {gene}\nDisease/Condition: {disease}\n"
        f"(Each evidence item has a [Ref N] tag — cite these in your report)\n"
    ]

    # カテゴリ別リファレンスレジストリ
    # category -> [(tag, citation_str), ...]
    ref_categories: dict[str, list[tuple[str, str]]] = {
        "paper":   [],   # PubMed 論文
        "disease": [],   # ClinVar / GWAS / OpenTargets
        "gene":    [],   # UniProt / IntAct / PubChem
        "drug":    [],   # ChEMBL
    }
    ref_counters: dict[str, int] = {k: 0 for k in ref_categories}

    CATEGORY_PREFIX = {
        "paper":   "Paper",
        "disease": None,   # ClinVar N / GWAS N / OpenTargets など個別プレフィックス
        "gene":    None,   # UniProt / IntAct など個別プレフィックス
        "drug":    "ChEMBL",
    }

    def add_ref(category: str, citation: str, prefix: str = "") -> str:
        """Register a citation and return its typed tag e.g. [Paper 1], [ClinVar 2]."""
        if not prefix:
            prefix = CATEGORY_PREFIX.get(category, "Ref")
        if prefix:
            ref_counters[category] += 1
            tag = f"[{prefix} {ref_counters[category]}]"
        else:
            tag = f"[{category}]"
        ref_categories[category].append((tag, citation))
        return tag

    # ── Gene/Protein info ──────────────────────────────────────────────────
    uni = ev.get("uniprot") or {}
    if uni and "error" not in uni:
        uniprot_id = uni.get("uniprot_id", gene)
        url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}" if uniprot_id else \
              "https://www.uniprot.org/"
        citation = (
            f"UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. "
            f"Nucleic Acids Res. 2023;51(D1):D523-D531. "
            f"Entry: {gene} ({uniprot_id}). {url}"
        )
        ref = add_ref("gene", citation, prefix="UniProt")
        sections.append(
            f"## Gene/Protein Information {ref}\n"
            f"- Protein: {uni.get('protein_name', 'N/A')}\n"
            f"- Function: {uni.get('function', 'N/A')[:500]}\n"
            f"- Subcellular location: {', '.join(uni.get('subcellular_location', [])) or 'N/A'}\n"
            f"- Keywords: {', '.join(uni.get('keywords', [])[:10]) or 'N/A'}\n"
            f"- GO terms (sample): {'; '.join(g['term'] for g in uni.get('go_terms', [])[:8]) or 'N/A'}\n"
        )

    # ── OpenTargets ────────────────────────────────────────────────────────
    ot = ev.get("opentargets") or {}
    if ot and "error" not in ot:
        ensg = ot.get("ensembl_id", "")
        efo  = ot.get("disease_id", "")
        ot_url = (
            f"https://platform.opentargets.org/evidence/{ensg}/{efo}"
            if ensg and efo else "https://platform.opentargets.org/"
        )
        citation = (
            f"Ochoa D, et al. Open Targets Platform: supporting systematic "
            f"drug-target identification and prioritisation. "
            f"Nucleic Acids Res. 2021;49(D1):D1302-D1310. "
            f"Target-Disease: {gene} ({ensg}) × {ot.get('disease_label', disease)} ({efo}). "
            f"{ot_url}"
        )
        ref = add_ref("disease", citation, prefix="OpenTargets")
        score = ot.get("association_score")
        score_str = f"{score:.3f}" if score is not None else "Not found"
        dt_scores = ot.get("datatype_scores", {})
        dt_str = "\n".join(f"  - {k}: {v:.3f}" for k, v in dt_scores.items()) or "  N/A"
        sections.append(
            f"## OpenTargets Association Evidence {ref}\n"
            f"- Overall association score: {score_str} (0–1)\n"
            f"- Evidence by data type:\n{dt_str}\n"
            f"- Disease label: {ot.get('disease_label', disease)}\n"
        )

    # ── GWAS ───────────────────────────────────────────────────────────────
    gwas_hits = ev.get("gwas") or []
    if gwas_hits:
        lines = []
        for h in gwas_hits[:5]:
            study_id = h.get("study_id", "")
            pub_date = h.get("pub_date", "")
            first_author = h.get("first_author", "")
            year = pub_date[:4] if pub_date else ""
            url = f"https://www.ebi.ac.uk/gwas/studies/{study_id}" if study_id else \
                  "https://www.ebi.ac.uk/gwas/"
            author_str = f"{first_author}, et al. " if first_author else ""
            year_str   = f"{year}. " if year else ""
            citation = (
                f"{author_str}GWAS Catalog Study {study_id}: {h.get('trait','')}. "
                f"{year_str}"
                f"Buniello A, et al. The NHGRI-EBI GWAS Catalog. "
                f"Nucleic Acids Res. 2019;47(D1):D1005-D1012. "
                f"{url}"
            )
            ref = add_ref("disease", citation, prefix="GWAS")
            snps = ", ".join(h.get("snps", [])[:2])
            lines.append(
                f"  - {h['trait']} | p={h['p_value']} | OR/Beta={h['or_beta']} "
                f"| SNPs: {snps} {ref}"
            )
        sections.append(f"## GWAS Associations (GWAS Catalog)\n" + "\n".join(lines) + "\n")

    # ── ClinVar ────────────────────────────────────────────────────────────
    cv_hits = ev.get("clinvar") or []
    if cv_hits:
        lines = []
        for v in cv_hits[:5]:
            var_id = str(v.get("uid", ""))
            url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{var_id}/" if var_id else \
                  f"https://www.ncbi.nlm.nih.gov/clinvar/?term={gene}"
            citation = (
                f"Landrum MJ, et al. ClinVar: improvements to accessing data. "
                f"Nucleic Acids Res. 2020;48(D1):D835-D844. "
                f"Variant: {v.get('title','')[:80]} "
                f"[Variation ID: {var_id}]. "
                f"Clinical significance: {v.get('clinical_significance','')}. "
                f"{url}"
            )
            ref = add_ref("disease", citation, prefix="ClinVar")
            lines.append(
                f"  - {v['title'][:70]} | {v['clinical_significance']} "
                f"| Condition: {v['condition']} {ref}"
            )
        sections.append(f"## ClinVar Pathogenic Variants\n" + "\n".join(lines) + "\n")

    # ── Existing drugs ─────────────────────────────────────────────────────
    drugs_chembl = ev.get("chembl") or []
    drugs_ot = ot.get("known_drugs", []) if ot else []
    all_drugs = {
        d.get("name") or d.get("drug", ""): d
        for d in drugs_chembl + drugs_ot
        if d.get("name") or d.get("drug")
    }
    if all_drugs:
        lines = []
        for name, d in list(all_drugs.items())[:8]:
            chembl_id = d.get("chembl_id", "")
            url = f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/" \
                  if chembl_id else "https://www.ebi.ac.uk/chembl/"
            phase = d.get("max_phase") or d.get("phase")
            mech  = d.get("mechanism", d.get("mechanism_of_action", ""))
            citation = (
                f"Mendez D, et al. ChEMBL: towards direct deposition of bioassay data. "
                f"Nucleic Acids Res. 2019;47(D1):D930-D940. "
                f"Compound: {name}"
                + (f" ({chembl_id})" if chembl_id else "")
                + (f". Max clinical phase: {phase}" if phase else "")
                + (f". Mechanism: {mech}" if mech else "")
                + f". {url}"
            )
            ref = add_ref("drug", citation, prefix="ChEMBL")
            lines.append(f"  - {name} | Phase: {phase} | Mechanism: {mech or 'N/A'} {ref}")
        sections.append(
            f"## Existing Drugs / Clinical Candidates Targeting {gene}\n"
            + "\n".join(lines) + "\n"
        )
    else:
        sections.append(
            f"## Existing Drugs\nNo approved drugs found targeting {gene} in ChEMBL/OpenTargets.\n"
        )

    # ── Protein interactions ───────────────────────────────────────────────
    interactions = ev.get("intact") or []
    if interactions:
        ia_url = f"https://www.ebi.ac.uk/intact/search?query={gene}"
        citation = (
            f"Orchard S, et al. The MIntAct project — IntAct as a common curation platform "
            f"for 11 molecular interaction databases. "
            f"Nucleic Acids Res. 2014;42(D1):D358-D363. "
            f"Query gene: {gene}. {ia_url}"
        )
        ref = add_ref("gene", citation, prefix="IntAct")
        partners = []
        for ix in interactions[:10]:
            partners.extend(ix.get("partners", []))
        unique_partners = list(dict.fromkeys(partners))[:10]
        sections.append(
            f"## Protein Interaction Network (IntAct) {ref}\n"
            f"- Key interactors: {', '.join(unique_partners) or 'N/A'}\n"
        )

    # ── Toxicity ───────────────────────────────────────────────────────────
    tox = ev.get("toxicity") or {}
    pb  = tox.get("pubchem_bioassay", {})
    ae  = tox.get("drug_adverse_events", {})
    ae_str = ""
    for drug_name, events in ae.items():
        ae_str += f"\n  {drug_name}: " + ", ".join(
            f"{e['reaction']}({e['count']})" for e in events[:3]
        )
    pc_url = f"https://pubchem.ncbi.nlm.nih.gov/#query={gene}&input_type=gene"
    citation_pc = (
        f"Kim S, et al. PubChem 2023 update. "
        f"Nucleic Acids Res. 2023;51(D1):D1373-D1380. "
        f"Gene bioassay query: {gene}. {pc_url}"
    )
    ref_pc = add_ref("gene", citation_pc, prefix="PubChem")
    sections.append(
        f"## Toxicity / Safety Signals\n"
        f"- PubChem BioAssay {ref_pc}: {pb.get('assay_count', 0)} assays\n"
        f"- Adverse events from related drugs:{ae_str or ' N/A'}\n"
        f"- Note: {tox.get('toxcast_note', '')}\n"
    )

    # ── Literature ─────────────────────────────────────────────────────────
    papers = ev.get("pubmed") or []
    if papers:
        paper_blocks = []
        for p in papers[:5]:
            pmid     = p.get("pmid", "")
            title    = p.get("title", "")
            journal  = p.get("journal", "")
            year     = p.get("year", "")
            authors  = p.get("authors", [])
            abstract = (p.get("abstract") or "").strip()
            url      = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""

            # Vancouver / NLM format citation
            if authors:
                author_str = ", ".join(authors[:3])
                if len(authors) > 3:
                    author_str += ", et al"
            else:
                author_str = "[Author(s) not available]"
            citation = (
                f"{author_str}. {title} "
                f"{journal}. {year}. "
                + (f"PMID: {pmid}. " if pmid else "")
                + url
            )
            ref = add_ref("paper", citation, prefix="Paper")

            # アブストラクトを 600 字で切り詰め（LLMへのコンテキスト量を制御）
            abstract_snippet = abstract[:600] + ("..." if len(abstract) > 600 else "") \
                               if abstract else "(Abstract not available)"

            paper_blocks.append(
                f"### {ref} [{year}] {title}\n"
                f"**Authors:** {author_str}  |  **Journal:** {journal}\n\n"
                f"**Abstract summary:**\n{abstract_snippet}\n"
            )
        sections.append(
            "## Recent Literature (PubMed) — with abstracts\n\n"
            + "\n".join(paper_blocks)
        )

    # ── gnomAD 制約スコア ──────────────────────────────────────────────────
    gnom = ev.get("gnomad") or {}
    if gnom and "error" not in gnom:
        url = gnom.get("url", "https://gnomad.broadinstitute.org/")
        citation = (
            f"Chen S, et al. A genomic mutational constraint map using variation in 76,156 human genomes. "
            f"Nature. 2024;625:92-100. "
            f"Gene: {gene} — pLI={gnom.get('pLI')}, LOEUF={gnom.get('LOEUF')}. {url}"
        )
        ref = add_ref("gene", citation, prefix="gnomAD")
        sections.append(
            f"## Population Genetics & Constraint (gnomAD) {ref}\n"
            f"- pLI: {gnom.get('pLI')} (>0.9 = highly constrained)\n"
            f"- LOEUF: {gnom.get('LOEUF')} (<0.35 = constrained)\n"
            f"- Missense O/E: {gnom.get('oe_missense')}\n"
            f"- Essentiality assessment: {gnom.get('essentiality')}\n"
        )

    # ── GTEx 組織発現 ──────────────────────────────────────────────────────
    gtex_data = ev.get("gtex") or {}
    if gtex_data and "error" not in gtex_data:
        url = gtex_data.get("url", "https://gtexportal.org/")
        citation = (
            f"GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. "
            f"Science. 2020;369(6509):1318-1330. "
            f"Gene expression: {gene}. {url}"
        )
        ref = add_ref("gene", citation, prefix="GTEx")
        top = gtex_data.get("top_tissues", [])[:5]
        top_str = ", ".join(f"{t['tissue']} ({t['tpm']:.1f} TPM)" for t in top)
        key = gtex_data.get("key_tissues", [])
        key_str = "\n".join(
            f"  - {t['tissue']}: {t['tpm']:.1f} TPM"
            for t in key if t["tpm"] > 0
        ) or "  N/A"
        sections.append(
            f"## Tissue Gene Expression (GTEx v8) {ref}\n"
            f"- Highest expression: {gtex_data.get('max_tissue')} ({gtex_data.get('max_tpm',0):.1f} TPM)\n"
            f"- Top 5 tissues: {top_str}\n"
            f"- Safety-relevant tissues:\n{key_str}\n"
        )

    # ── Human Protein Atlas ────────────────────────────────────────────────
    hpa_data = ev.get("hpa") or {}
    if hpa_data and "error" not in hpa_data:
        url = hpa_data.get("url", "https://www.proteinatlas.org/")
        citation = (
            f"Uhlén M, et al. Tissue-based map of the human proteome. "
            f"Science. 2015;347(6220):1260419. "
            f"Gene: {gene}. {url}"
        )
        ref = add_ref("gene", citation, prefix="HPA")
        subcell = ", ".join(hpa_data.get("subcellular", [])[:5]) or "N/A"
        prot_class = ", ".join(hpa_data.get("protein_class", [])[:5]) or "N/A"
        tissue_lines = "\n".join(
            f"  - {t['tissue']} / {t['cell_type']}: {t['level']} ({t['reliability']})"
            for t in hpa_data.get("tissue_expression", [])[:8]
        ) or "  N/A"
        sections.append(
            f"## Protein Expression Profile (Human Protein Atlas) {ref}\n"
            f"- Protein class: {prot_class}\n"
            f"- Subcellular location: {subcell}\n"
            f"- High/Medium expression tissues:\n{tissue_lines}\n"
        )

    # ── DGIdb 薬剤-遺伝子相互作用 ─────────────────────────────────────────
    dgi_data = ev.get("dgidb") or []
    if dgi_data:
        url = f"https://dgidb.org/genes/{gene}#interactions"
        citation = (
            f"Cannon M, et al. DGIdb 4.0: updates to the drug-gene interaction database. "
            f"Nucleic Acids Res. 2024;52(D1):D1227-D1235. "
            f"Gene: {gene}. {url}"
        )
        ref = add_ref("drug", citation, prefix="DGIdb")
        approved = [d for d in dgi_data if d.get("approved")]
        lines = []
        for d in dgi_data[:8]:
            app_str = "✓ Approved" if d.get("approved") else "Investigational"
            i_type  = d.get("interaction_type", "")
            direc   = d.get("directionality", "")
            lines.append(
                f"  - {d['drug_name']} | {app_str} | {i_type} {direc} "
                f"| Sources: {', '.join(d.get('sources',[])[:3])}"
            )
        sections.append(
            f"## Drug-Gene Interactions (DGIdb) {ref}\n"
            f"- Total: {len(dgi_data)} interactions ({len(approved)} approved drugs)\n"
            + "\n".join(lines) + "\n"
        )

    # ── ClinicalTrials.gov ─────────────────────────────────────────────────
    ct_data = ev.get("clinicaltrials") or []
    if ct_data:
        citation = (
            f"ClinicalTrials.gov. U.S. National Library of Medicine. "
            f"Search: {gene} × {disease}. "
            f"https://clinicaltrials.gov/search?cond={disease.replace(' ','%20')}&intr={gene}"
        )
        ref = add_ref("disease", citation, prefix="ClinicalTrials")
        lines = []
        for t in ct_data[:6]:
            ivs = ", ".join(t.get("interventions", [])[:3])
            lines.append(
                f"  - {t['nct_id']} | {t['phase']} | {t['status']} | {t['title'][:60]} | Interventions: {ivs}"
            )
        sections.append(
            f"## Clinical Trials (ClinicalTrials.gov) {ref}\n"
            f"- {len(ct_data)} trials found for {gene} × {disease}\n"
            + "\n".join(lines) + "\n"
        )

    # ── AlphaFold 構造情報 ─────────────────────────────────────────────────
    af_data = ev.get("alphafold") or {}
    if af_data and "error" not in af_data:
        url = af_data.get("view_url", "https://alphafold.ebi.ac.uk/")
        citation = (
            f"Jumper J, et al. Highly accurate protein structure prediction with AlphaFold. "
            f"Nature. 2021;596:583-589. "
            f"Entry: {af_data.get('entry_id','')} (UniProt: {af_data.get('uniprot_id','')}). {url}"
        )
        ref = add_ref("gene", citation, prefix="AlphaFold")
        sections.append(
            f"## Predicted 3D Structure (AlphaFold) {ref}\n"
            f"- Entry: {af_data.get('entry_id','')}\n"
            f"- Mean pLDDT: {af_data.get('mean_plddt')} / 100\n"
            f"- Structural confidence: {af_data.get('confidence')}\n"
            f"- Structure URL: {af_data.get('pdb_url','')}\n"
        )

    # ── Reactome パスウェイ ────────────────────────────────────────────────
    react_data = ev.get("reactome") or []
    if react_data:
        url = f"https://reactome.org/content/query?q={gene}&species=Homo+sapiens"
        citation = (
            f"Milacic M, et al. The Reactome Pathway Knowledgebase 2024. "
            f"Nucleic Acids Res. 2024;52(D1):D672-D678. "
            f"Gene: {gene}. {url}"
        )
        ref = add_ref("gene", citation, prefix="Reactome")
        disease_pw = [p for p in react_data if p.get("is_disease")]
        all_pw     = react_data[:10]
        lines = [f"  - {p['name']} [{p['pathway_id']}]" for p in all_pw]
        sections.append(
            f"## Reactome Pathways {ref}\n"
            f"- Total pathways: {len(react_data)} "
            f"({len(disease_pw)} disease-annotated)\n"
            + "\n".join(lines) + "\n"
        )

    # ── Reference list（カテゴリ別） ──────────────────────────────────────
    ref_section_lines = ["## References",
                         "(Cite inline using the tags below. e.g. [Paper 1], [ClinVar 2], [OpenTargets 1])\n"]

    CATEGORY_HEADERS = {
        "paper":   "### Papers (PubMed)",
        "disease": "### Disease & Genetic Databases  (ClinVar / GWAS / OpenTargets / ClinicalTrials)",
        "gene":    "### Gene & Protein Databases  (UniProt / IntAct / gnomAD / GTEx / HPA / AlphaFold / Reactome / PubChem)",
        "drug":    "### Drug Databases  (ChEMBL / DGIdb)",
    }
    any_ref = False
    for cat, header in CATEGORY_HEADERS.items():
        entries = ref_categories.get(cat, [])
        if entries:
            any_ref = True
            ref_section_lines.append(header)
            for tag, cite in entries:
                ref_section_lines.append(f"{tag} {cite}")
            ref_section_lines.append("")

    if any_ref:
        sections.append("\n".join(ref_section_lines))

    return "\n".join(sections)
