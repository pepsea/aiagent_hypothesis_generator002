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
    """Convert aggregated evidence into a structured text context for LLM."""
    gene = aggregated["gene"]
    disease = aggregated["disease"]
    ev = aggregated["evidence"]

    sections = [f"# Evidence Summary for Drug Target Hypothesis\n\nGene: {gene}\nDisease/Condition: {disease}\n"]

    # --- Gene/Protein info ---
    uni = ev.get("uniprot") or {}
    if uni and "error" not in uni:
        sections.append(f"""## Gene/Protein Information
- Protein: {uni.get('protein_name', 'N/A')}
- Function: {uni.get('function', 'N/A')[:500]}
- Subcellular location: {', '.join(uni.get('subcellular_location', [])) or 'N/A'}
- Keywords: {', '.join(uni.get('keywords', [])[:10]) or 'N/A'}
- GO terms (sample): {'; '.join(f"{g['term']}" for g in uni.get('go_terms', [])[:8]) or 'N/A'}
""")

    # --- OpenTargets evidence ---
    ot = ev.get("opentargets") or {}
    if ot and "error" not in ot:
        score = ot.get("association_score")
        score_str = f"{score:.3f}" if score is not None else "Not found"
        dt_scores = ot.get("datatype_scores", {})
        dt_str = "\n".join(f"  - {k}: {v:.3f}" for k, v in dt_scores.items()) if dt_scores else "  N/A"

        sections.append(f"""## OpenTargets Association Evidence
- Overall association score: {score_str} (scale 0–1, higher = stronger)
- Evidence by data type:
{dt_str}
- Disease label in OpenTargets: {ot.get('disease_label', disease)}
""")

    # --- GWAS ---
    gwas_hits = ev.get("gwas") or []
    if gwas_hits:
        gwas_str = "\n".join(
            f"  - {h['trait']} | p={h['p_value']} | OR/Beta={h['or_beta']} | SNPs: {', '.join(h['snps'][:2])}"
            for h in gwas_hits[:5]
        )
        sections.append(f"""## GWAS Associations (GWAS Catalog)
{gwas_str}
""")

    # --- ClinVar ---
    cv_hits = ev.get("clinvar") or []
    if cv_hits:
        cv_str = "\n".join(
            f"  - {v['title']} | {v['clinical_significance']} | Condition: {v['condition']}"
            for v in cv_hits[:5]
        )
        sections.append(f"""## ClinVar Pathogenic Variants
{cv_str}
""")

    # --- Existing drugs ---
    drugs_chembl = ev.get("chembl") or []
    drugs_ot = ot.get("known_drugs", []) if ot else []

    all_drugs = {d.get("name") or d.get("drug", ""): d for d in drugs_chembl + drugs_ot if d.get("name") or d.get("drug")}
    if all_drugs:
        drug_str = "\n".join(
            f"  - {name} | Phase: {d.get('max_phase') or d.get('phase')} | Mechanism: {d.get('mechanism', d.get('mechanism_of_action', 'N/A'))}"
            for name, d in list(all_drugs.items())[:8]
        )
        sections.append(f"""## Existing Drugs / Clinical Candidates Targeting {gene}
{drug_str}
""")
    else:
        sections.append(f"## Existing Drugs\nNo approved drugs found targeting {gene} in ChEMBL/OpenTargets.\n")

    # --- Protein interactions ---
    interactions = ev.get("intact") or []
    if interactions:
        partners = []
        for ix in interactions[:10]:
            partners.extend(ix.get("partners", []))
        unique_partners = list(dict.fromkeys(partners))[:10]
        sections.append(f"""## Protein Interaction Network (IntAct)
- Key interactors: {', '.join(unique_partners) or 'N/A'}
""")

    # --- Toxicity ---
    tox = ev.get("toxicity") or {}
    pb = tox.get("pubchem_bioassay", {})
    ae = tox.get("drug_adverse_events", {})
    ae_str = ""
    for drug_name, events in ae.items():
        ae_str += f"\n  {drug_name}: " + ", ".join(f"{e['reaction']}({e['count']})" for e in events[:3])

    sections.append(f"""## Toxicity / Safety Signals
- PubChem BioAssay: {pb.get('assay_count', 0)} assays involving {gene}
- Adverse events from related drugs:{ae_str or ' N/A'}
- Note: {tox.get('toxcast_note', '')}
""")

    # --- Literature ---
    papers = ev.get("pubmed") or []
    if papers:
        paper_str = "\n".join(
            f"  - [{p['year']}] {p['title']} ({p['journal']})"
            for p in papers[:5]
        )
        sections.append(f"""## Recent Literature (PubMed)
{paper_str}
""")

    return "\n".join(sections)
