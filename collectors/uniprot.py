"""UniProt REST API (CC BY 4.0) — gene/protein info."""
import requests

BASE = "https://rest.uniprot.org/uniprotkb"

def get_protein_info(gene_symbol: str, organism: str = "human") -> dict:
    """Return key protein annotations for a human gene."""
    query = f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true"
    r = requests.get(f"{BASE}/search", params={
        "query": query,
        "fields": "accession,gene_names,protein_name,organism_name,cc_function,"
                  "cc_subcellular_location,cc_disease,go,keyword,length,ft_domain",
        "format": "json",
        "size": 1,
    }, timeout=15)
    r.raise_for_status()
    results = r.json().get("results", [])
    if not results:
        return {"error": f"No UniProt entry for {gene_symbol}"}

    entry = results[0]
    acc = entry.get("primaryAccession", "")

    # Extract fields
    function_cc = ""
    subcell = []
    diseases = []
    for comment in entry.get("comments", []):
        ctype = comment.get("commentType", "")
        if ctype == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                function_cc = texts[0].get("value", "")
        elif ctype == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                subcell.append(loc.get("location", {}).get("value", ""))
        elif ctype == "DISEASE":
            d = comment.get("disease", {})
            diseases.append({
                "name": d.get("diseaseId", ""),
                "description": d.get("description", ""),
            })

    go_terms = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "GO":
            go_terms.append({
                "id": xref.get("id", ""),
                "term": next((p.get("value") for p in xref.get("properties", [])
                              if p.get("key") == "GoTerm"), ""),
                "aspect": next((p.get("value") for p in xref.get("properties", [])
                                if p.get("key") == "GoEvidenceType"), ""),
            })

    keywords = [kw.get("name", "") for kw in entry.get("keywords", [])]

    protein_name = ""
    pn = entry.get("proteinDescription", {})
    if pn.get("recommendedName"):
        protein_name = pn["recommendedName"].get("fullName", {}).get("value", "")

    return {
        "accession": acc,
        "gene": gene_symbol,
        "protein_name": protein_name,
        "function": function_cc,
        "subcellular_location": subcell,
        "associated_diseases": diseases,
        "go_terms": go_terms[:15],
        "keywords": keywords[:20],
        "sequence_length": entry.get("sequence", {}).get("length"),
    }


def _first_function_text(entry: dict) -> str:
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                return texts[0].get("value", "")
    return ""


def _gene_symbol_of(entry: dict) -> str:
    for g in entry.get("genes", []):
        name = (g.get("geneName") or {}).get("value")
        if name:
            return name
    return ""


def get_functions_for_genes(gene_symbols: list[str]) -> dict[str, dict]:
    """複数遺伝子の機能情報を1回のバッチクエリで取得する（PPIパートナー用）。

    Returns:
        {GENE_SYMBOL(upper): {"protein_name", "function", "accession"}}
    """
    genes = [g for g in dict.fromkeys(g.strip() for g in gene_symbols) if g]
    if not genes:
        return {}

    # (gene_exact:A OR gene_exact:B ...) AND ヒト・reviewed のバッチ検索
    or_clause = " OR ".join(f"gene_exact:{g}" for g in genes)
    query = f"({or_clause}) AND organism_id:9606 AND reviewed:true"
    try:
        r = requests.get(f"{BASE}/search", params={
            "query": query,
            "fields": "accession,gene_names,protein_name,cc_function",
            "format": "json",
            "size": min(len(genes) * 2, 500),
        }, timeout=25)
        r.raise_for_status()
        results = r.json().get("results", [])
    except Exception:
        return {}

    out: dict[str, dict] = {}
    for entry in results:
        sym = _gene_symbol_of(entry).upper()
        if not sym or sym in out:
            continue
        pn = entry.get("proteinDescription", {})
        protein_name = ""
        if pn.get("recommendedName"):
            protein_name = pn["recommendedName"].get("fullName", {}).get("value", "")
        out[sym] = {
            "protein_name": protein_name,
            "function": _first_function_text(entry),
            "accession": entry.get("primaryAccession", ""),
        }
    return out
