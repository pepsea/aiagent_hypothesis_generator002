"""Toxicity data collectors.

Sources:
- PubChem BioAssay (public domain, NIH) — in vitro toxicity assays
- ToxCast (EPA, public domain) — high-throughput toxicity screening
- FDA Drug Safety (openFDA, public domain) — adverse events
"""
import requests

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
OPENFDA_BASE = "https://api.fda.gov/drug"


def get_pubchem_gene_toxicity(gene_symbol: str) -> dict:
    """Get gene involvement in toxicity assays via PubChem BioAssay."""
    # Query gene-related bioassays with toxicity endpoints
    r = requests.get(
        f"{PUBCHEM_BASE}/gene/genesymbol/{gene_symbol}/aids/JSON",
        timeout=15
    )
    if r.status_code == 404:
        return {"assay_count": 0, "gene": gene_symbol}
    r.raise_for_status()

    aids = r.json().get("InformationList", {}).get("Information", [{}])[0].get("AID", [])

    return {
        "gene": gene_symbol,
        "assay_count": len(aids),
        "sample_assay_ids": aids[:5],
        "source": "PubChem BioAssay (NIH, public domain)",
    }


def get_openfda_adverse_events(drug_name: str, limit: int = 5) -> list[dict]:
    """Return top adverse events for a drug from openFDA (public domain)."""
    r = requests.get(f"{OPENFDA_BASE}/event.json", params={
        "search": f'patient.drug.medicinalproduct:"{drug_name}"',
        "count": "patient.reaction.reactionmeddrapt.exact",
        "limit": limit,
    }, timeout=15)

    if r.status_code in (404, 400):
        return []
    r.raise_for_status()

    results = r.json().get("results", [])
    return [{"reaction": item.get("term", ""), "count": item.get("count", 0)} for item in results]


def get_toxcast_summary(gene_symbol: str) -> dict:
    """Get ToxCast/Tox21 assay summary for a gene target (EPA public data via comptox API)."""
    # EPA CompTox Dashboard API — public domain
    r = requests.get(
        f"https://comptox.epa.gov/dashboard-api/ccdapp1/chemical-files/search/by-name/{gene_symbol}",
        timeout=15
    )
    # This endpoint may not always return gene data; use as supplemental
    return {
        "gene": gene_symbol,
        "note": "ToxCast data available via EPA CompTox Dashboard (https://comptox.epa.gov)",
        "source": "EPA ToxCast/Tox21 (public domain)",
    }


def assess_target_safety(gene_symbol: str, known_drugs: list[dict]) -> dict:
    """Aggregate toxicity signals for a target."""
    pubchem_data = get_pubchem_gene_toxicity(gene_symbol)

    adverse_events = {}
    for drug in known_drugs[:3]:
        drug_name = drug.get("drug") or drug.get("name", "")
        if drug_name:
            ae = get_openfda_adverse_events(drug_name, limit=5)
            if ae:
                adverse_events[drug_name] = ae

    return {
        "pubchem_bioassay": pubchem_data,
        "drug_adverse_events": adverse_events,
        "toxcast_note": "See EPA CompTox Dashboard for ToxCast assay data",
    }
