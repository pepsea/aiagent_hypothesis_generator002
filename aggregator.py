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

    # Numbered reference registry: ref_id -> {source, id, url, label}
    refs: list[dict] = []

    def add_ref(source: str, label: str, url: str = "", ref_id: str = "") -> str:
        """Register a reference and return its [Ref N] tag."""
        n = len(refs) + 1
        refs.append({"n": n, "source": source, "id": ref_id, "url": url, "label": label})
        return f"[Ref {n}]"

    # ── Gene/Protein info ──────────────────────────────────────────────────
    uni = ev.get("uniprot") or {}
    if uni and "error" not in uni:
        uniprot_id = uni.get("uniprot_id", gene)
        uniprot_url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}" if uniprot_id else ""
        ref = add_ref("UniProt", f"{gene} protein entry", uniprot_url, uniprot_id)
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
            if ensg and efo else ""
        )
        ref = add_ref("OpenTargets", f"{gene}×{disease} association", ot_url, f"{ensg}/{efo}")
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
            url = f"https://www.ebi.ac.uk/gwas/studies/{study_id}" if study_id else \
                  "https://www.ebi.ac.uk/gwas/"
            ref = add_ref("GWAS Catalog", h.get("trait", ""), url, study_id)
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
            ref = add_ref("ClinVar", v.get("title", "")[:60], url, var_id)
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
            ref = add_ref("ChEMBL", name, url, chembl_id)
            phase = d.get("max_phase") or d.get("phase")
            mech  = d.get("mechanism", d.get("mechanism_of_action", "N/A"))
            lines.append(f"  - {name} | Phase: {phase} | Mechanism: {mech} {ref}")
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
        ref = add_ref("IntAct", f"{gene} interaction network", ia_url, gene)
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
    ref_pc = add_ref("PubChem BioAssay", f"{gene} bioassay data", pc_url, gene)
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
            pmid = p.get("pmid", "")
            url  = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
            ref  = add_ref("PubMed", p.get("title", "")[:80], url, pmid)
            lines.append(
                f"  - [{p['year']}] {p['title'][:80]} ({p['journal'][:40]}) {ref}"
            )
        sections.append(f"## Recent Literature (PubMed)\n" + "\n".join(lines) + "\n")

    # ── Reference list ─────────────────────────────────────────────────────
    if refs:
        ref_lines = []
        for r in refs:
            url_part = f" — {r['url']}" if r["url"] else ""
            id_part  = f" (ID: {r['id']})" if r["id"] else ""
            ref_lines.append(f"[Ref {r['n']}] {r['source']}: {r['label']}{id_part}{url_part}")
        sections.append(
            "## References\n" + "\n".join(ref_lines) + "\n\n"
            "(When writing the report, cite evidence using [Ref N] inline.)\n"
        )

    return "\n".join(sections)
