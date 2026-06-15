"""Aggregate multi-source evidence into a structured context for LLM."""
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from collectors import pubmed, opentargets, intact, uniprot, gwas, chembl, toxicity


def collect_all(
    gene: str,
    disease: str,
    verbose: bool = True,
    gene_id: str = None,
    disease_id: str = None,
) -> dict:
    """Run all collectors in parallel and return aggregated evidence.

    gene_id / disease_id: optional Ensembl / EFO IDs from OpenTargets widget.
    When provided, they skip redundant ID-resolution queries inside opentargets.py.
    """

    def log(msg):
        if verbose:
            print(f"  [+] {msg}")

    tasks = {
        "pubmed": lambda: pubmed.search_pubmed(gene, disease, max_results=5),
        "opentargets": lambda: opentargets.get_target_disease_evidence(
            gene, disease, gene_id=gene_id, disease_id=disease_id
        ),
        "uniprot": lambda: uniprot.get_protein_info(gene),
        "intact": lambda: intact.get_interactions(gene, max_results=15),
        "gwas": lambda: gwas.get_gwas_associations(gene, disease),
        "clinvar": lambda: gwas.get_clinvar_variants(gene),
        "chembl": lambda: chembl.get_drugs_for_target(gene),
    }

    results = {}
    errors = {}

    with ThreadPoolExecutor(max_workers=6) as executor:
        future_to_key = {executor.submit(fn): key for key, fn in tasks.items()}
        for future in as_completed(future_to_key):
            key = future_to_key[future]
            try:
                results[key] = future.result()
                log(f"{key}: OK")
            except Exception as e:
                errors[key] = str(e)
                results[key] = None
                log(f"{key}: FAILED ({e})")

    # Toxicity requires known drugs from chembl/opentargets
    known_drugs = []
    if results.get("chembl"):
        known_drugs.extend(results["chembl"])
    if results.get("opentargets") and isinstance(results["opentargets"], dict):
        known_drugs.extend(results["opentargets"].get("known_drugs", []))

    try:
        results["toxicity"] = toxicity.assess_target_safety(gene, known_drugs)
        log("toxicity: OK")
    except Exception as e:
        errors["toxicity"] = str(e)
        results["toxicity"] = None
        log(f"toxicity: FAILED ({e})")

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

    # Numbered reference registry; each entry holds a fully-formatted citation string.
    refs: list[str] = []

    def add_ref(citation: str) -> str:
        """Register a formatted citation and return its [Ref N] tag."""
        n = len(refs) + 1
        refs.append(citation)
        return f"[Ref {n}]"

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
        ref = add_ref(citation)
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
        ref = add_ref(citation)
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
            ref = add_ref(citation)
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
            ref = add_ref(citation)
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
            ref = add_ref(citation)
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
        ref = add_ref(citation)
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
    ref_pc = add_ref(citation_pc)
    sections.append(
        f"## Toxicity / Safety Signals\n"
        f"- PubChem BioAssay {ref_pc}: {pb.get('assay_count', 0)} assays\n"
        f"- Adverse events from related drugs:{ae_str or ' N/A'}\n"
        f"- Note: {tox.get('toxcast_note', '')}\n"
    )

    # ── Literature ─────────────────────────────────────────────────────────
    papers = ev.get("pubmed") or []
    if papers:
        lines = []
        for p in papers[:5]:
            pmid    = p.get("pmid", "")
            title   = p.get("title", "")
            journal = p.get("journal", "")
            year    = p.get("year", "")
            authors = p.get("authors", [])
            url     = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""

            # Vancouver / NLM format
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
            ref = add_ref(citation)
            lines.append(
                f"  - [{year}] {title[:80]} ({journal[:40]}) {ref}"
            )
        sections.append(f"## Recent Literature (PubMed)\n" + "\n".join(lines) + "\n")

    # ── Reference list ─────────────────────────────────────────────────────
    if refs:
        ref_lines = [f"[Ref {i+1}] {cite}" for i, cite in enumerate(refs)]
        sections.append(
            "## References\n" + "\n".join(ref_lines) + "\n\n"
            "(When writing the report, cite evidence using [Ref N] inline.)\n"
        )

    return "\n".join(sections)
