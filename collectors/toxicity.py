"""Toxicity data collectors.

Sources:
- ToxCast (EPA CTX Bioactivity API, public domain) — high-throughput toxicity
  screening, gene-level assay activity. License: free for non-commercial and
  commercial use. Requires a personal API key (request via ccte_api@epa.gov),
  set as the CCTE_API_KEY environment variable.
- FDA Drug Safety (openFDA, public domain) — adverse events
"""
import os

import requests

OPENFDA_BASE = "https://api.fda.gov/drug"
CTX_BASE = "https://comptox.epa.gov/ctx-api"


def get_toxcast_gene_assays(gene_symbol: str) -> dict:
    """Get ToxCast assay activity summary for a gene target via the EPA CTX Bioactivity API.

    Requires CCTE_API_KEY env var (request a free key at ccte_api@epa.gov).
    """
    api_key = os.environ.get("CCTE_API_KEY")
    if not api_key:
        return {
            "gene": gene_symbol,
            "available": False,
            "note": "CCTE_API_KEY 未設定のため ToxCast データ取得不可 "
                    "(ccte_api@epa.gov にキーを申請してください)",
        }

    r = requests.get(
        f"{CTX_BASE}/bioactivity/assay/search/by-gene/{gene_symbol}",
        headers={"x-api-key": api_key},
        timeout=15,
    )
    if r.status_code == 404:
        return {"gene": gene_symbol, "available": True, "assay_count": 0, "assays": []}
    r.raise_for_status()

    data = r.json()
    assays = data if isinstance(data, list) else [data]
    assays = [a for a in assays if a]

    return {
        "gene": gene_symbol,
        "available": True,
        "assay_count": len(assays),
        "assays": assays[:10],
        "source": "EPA ToxCast/Tox21 via CTX Bioactivity API (public domain)",
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


def assess_target_safety(gene_symbol: str, known_drugs: list[dict]) -> dict:
    """Aggregate toxicity signals for a target."""
    toxcast_data = get_toxcast_gene_assays(gene_symbol)

    adverse_events = {}
    for drug in known_drugs[:3]:
        drug_name = drug.get("drug") or drug.get("name", "")
        if drug_name:
            ae = get_openfda_adverse_events(drug_name, limit=5)
            if ae:
                adverse_events[drug_name] = ae

    return {
        "toxcast": toxcast_data,
        "drug_adverse_events": adverse_events,
    }
