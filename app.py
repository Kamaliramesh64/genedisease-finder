"""GeneDisease Finder - single-file Flask app with all four modules.

    Module 1  Input, validation, gene features (structural and functional)
    Module 2  Pathways, common features, interaction network
    Module 3  Disease ranking, evidence filter, related diseases, drug targets
    Module 4  Compare two gene sets, charts, network views, CSV and PDF reports

Run:
    pip install flask pandas openpyxl requests
    python app.py
    open http://127.0.0.1:5000

Data sources (all free, no login):
    MyGene.info, UniProt, AlphaFold DB   gene and protein features
    Enrichr (Reactome library)           pathway enrichment
    STRING                               protein interaction network
    Open Targets Platform                gene-disease associations, drug targets, tissue expression
    MyVariant.info (ClinVar)             known variants per gene
"""
import csv
import io
import math
import os
import re
import uuid
from collections import OrderedDict, defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache

import pandas as pd
import requests
from flask import Flask, Response, abort, flash, redirect, render_template_string, request, url_for

# ----------------------------------------------------------------- settings

MYGENE_URL = "https://mygene.info/v3"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"
ALPHAFOLD_ENTRY_URL = "https://alphafold.ebi.ac.uk/entry"
ENRICHR_URL = "https://maayanlab.cloud/Enrichr"
ENRICHR_LIBRARIES = ["Reactome_Pathways_2024", "Reactome_2022"]  # first one that returns rows is used
STRING_URL = "https://string-db.org/api/json/network"
OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
MYVARIANT_URL = "https://myvariant.info/v1/query"

TIMEOUT = 30
BATCH_SIZE = 500
MAX_GENES = 5000            # most genes accepted in one list
MAX_FEATURE_GENES = 30      # gene feature cards shown
MAX_DISEASE_GENES = 100     # genes sent to Open Targets (one request each)
MAX_NETWORK_GENES = 100     # genes sent to STRING
TOP_PATHWAYS = 15
TOP_DISEASES = 15
MAX_GO_TERMS = 8
SAVED_RESULTS = 20          # results kept in memory for the download buttons

HEADER_NAMES = {"gene", "genes", "symbol", "symbols", "gene_symbol", "gene symbol", "gene_name", "gene name"}
UNIPROT_FEATURE_TYPES = {"Domain", "DNA binding", "Zinc finger", "Motif"}
GENE_FIELDS = "symbol,name,summary,type_of_gene,map_location,genomic_pos,exons,uniprot,go,ensembl.gene"
DRUGGABLE_LABELS = {
    "Approved Drug", "Advanced Clinical", "Phase 1 Clinical", "Structure with Ligand",
    "High-Quality Ligand", "High-Quality Pocket", "Med-Quality Pocket", "Druggable Family",
}
GENERIC_DISEASES = {          # hidden when "Hide very generic diseases" is ticked
    "cancer", "neoplasm", "carcinoma", "tumor", "tumour", "malignant neoplasm", "benign neoplasm",
    "disease", "genetic disease", "syndrome", "carcinoma in situ", "malignant tumor",
}
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]
SAMPLE_GENES = "TP53 BRCA1 EGFR KRAS PTEN PIK3CA AKT1 MYC CDKN2A ERBB2"
SAMPLE_GENES_B = "IL6 TNF SERPINE1 INS PPARG AKT1 PTEN"

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024  # 5 MB upload limit
app.secret_key = os.environ.get("SECRET_KEY", "change-this-before-deploying")
RESULTS = OrderedDict()


class ApiError(Exception):
    """An online database could not be reached or returned something unusable."""


def http(method, url, **kwargs):
    """Single place for every web request, so errors and timeouts are handled once."""
    try:
        response = requests.request(method, url, timeout=TIMEOUT, **kwargs)
        response.raise_for_status()
        return response.json()
    except (requests.RequestException, ValueError) as exc:
        detail = str(exc)
        reply = getattr(exc, "response", None)
        if reply is not None:
            detail = f"HTTP {reply.status_code} {reply.text[:150]}"
        raise ApiError(f"{url.split('/')[2]}: {detail}") from exc


# =========================================================== MODULE 1: input

def clean_symbols(items):
    """Uppercase, strip, drop blanks and duplicates, keep the original order."""
    seen, out = set(), []
    for item in items:
        if item is None:
            continue
        symbol = str(item).strip().strip("\"'").upper()
        if not symbol or symbol in {"NAN", "NONE"} or symbol in seen:
            continue
        seen.add(symbol)
        out.append(symbol)
    return out


def parse_genes_text(text):
    """Split pasted text on commas, semicolons, spaces and new lines."""
    return clean_symbols(re.split(r"[,;\s]+", text)) if text else []


def _column_from_table(df):
    """Use the column named like 'gene' or 'symbol' if there is one, else the first column."""
    if df.empty:
        return []
    first_row = [str(v).strip().lower() for v in df.iloc[0]]
    for idx, value in enumerate(first_row):
        if value in HEADER_NAMES:
            return list(df.iloc[1:, idx])
    return list(df.iloc[:, 0])


def _delimiter(text, filename):
    if filename.endswith(".tsv"):
        return "\t"
    try:
        return csv.Sniffer().sniff(text[:2048], delimiters=",;\t|").delimiter
    except csv.Error:
        return ","


def parse_genes_file(file_storage):
    """Read genes from an uploaded .csv, .tsv, .txt, .xlsx or .xls file."""
    if file_storage is None or not getattr(file_storage, "filename", ""):
        return []
    name = file_storage.filename.lower()
    data = file_storage.read()
    if name.endswith((".xlsx", ".xls")):
        df = pd.read_excel(io.BytesIO(data), header=None, dtype=str)
        return clean_symbols(_column_from_table(df))
    text = data.decode("utf-8-sig", errors="ignore")
    if name.endswith((".csv", ".tsv")):
        df = pd.read_csv(io.StringIO(text), header=None, dtype=str, sep=_delimiter(text, name), engine="python")
        return clean_symbols(_column_from_table(df))
    if name.endswith(".txt"):
        return parse_genes_text(text)
    raise ValueError("Unsupported file type. Use .csv, .tsv, .txt or .xlsx")


# ----------------------------------------------------------- validation

def suggest_symbol(symbol):
    """Closest official human symbol that starts with the typed text, or None."""
    try:
        result = http("GET", f"{MYGENE_URL}/query",
                      params={"q": f"symbol:{symbol}*", "species": "human", "fields": "symbol", "size": 1})
    except ApiError:
        return None
    hits = result.get("hits", [])
    return hits[0].get("symbol") if hits else None


def _ensembl_id(hit):
    ensembl = hit.get("ensembl")
    if isinstance(ensembl, list):
        ensembl = ensembl[0] if ensembl else None
    return ensembl.get("gene") if isinstance(ensembl, dict) else None


def _chromosome(hit):
    """Chromosome name (1 to 22, X, Y) from a MyGene.info hit, or None."""
    pos = hit.get("genomic_pos")
    for item in (pos if isinstance(pos, list) else [pos]):
        chrom = str((item or {}).get("chr", "")).upper()
        if chrom in CHROMOSOMES:
            return chrom
    return None


def _lookup(batch, scopes):
    """MyGene.info lookup of symbols in one scope; returns {typed symbol: hit} for the genes found."""
    hits = http("POST", f"{MYGENE_URL}/query", data={
        "q": ",".join(batch), "scopes": scopes, "species": "human",
        "fields": "symbol,entrezgene,ensembl.gene,genomic_pos.chr", "size": 1,
    })
    found = {}
    for hit in hits:
        query = str(hit.get("query", "")).upper()
        if hit.get("notfound") or "entrezgene" not in hit or query in found:
            continue
        found[query] = hit
    return found


def validate_symbols(symbols):
    """Check symbols against MyGene.info (human genes).

    Returns (valid, invalid, notes):
        valid   {typed symbol: {"symbol": official, "entrez": id, "ensembl": id}}
        invalid [{"symbol": typed, "suggestion": official or None}]
        notes   messages such as 'HER2 was matched to the official symbol ERBB2.'
    """
    valid, invalid, notes = {}, [], []
    for start in range(0, len(symbols), BATCH_SIZE):
        batch = symbols[start:start + BATCH_SIZE]
        found = _lookup(batch, "symbol")                                   # official symbols first
        rest = [sym for sym in batch if sym not in found]
        if rest:
            found.update(_lookup(rest, "alias"))                           # then aliases for what is left
        for query, hit in found.items():
            official = str(hit.get("symbol", query)).upper()
            valid[query] = {"symbol": official, "entrez": str(hit["entrezgene"]), "ensembl": _ensembl_id(hit), "chr": _chromosome(hit)}
            if official != query:
                notes.append(f"{query} was matched to the official symbol {official}.")
        for symbol in batch:
            if symbol not in found:
                invalid.append({"symbol": symbol, "suggestion": suggest_symbol(symbol)})
    return valid, invalid, notes


# ----------------------------------------------------- gene features

def _as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _first_swissprot(uniprot_field):
    if not isinstance(uniprot_field, dict):
        return None
    accession = _as_list(uniprot_field.get("Swiss-Prot"))
    return accession[0] if accession else None


def _go_terms(go_field, branch):
    terms, seen = [], set()
    for item in _as_list((go_field or {}).get(branch)):
        term = item.get("term") if isinstance(item, dict) else None
        if term and term not in seen:
            seen.add(term)
            terms.append(term)
    return terms[:MAX_GO_TERMS]


def _exon_summary(exons_field):
    transcripts = _as_list(exons_field)
    counts = [len(t.get("position", [])) for t in transcripts if isinstance(t, dict)]
    return len(transcripts), (max(counts) if counts else None)


def parse_uniprot(data):
    length = (data.get("sequence") or {}).get("length")
    domains = []
    for feature in data.get("features", []):
        if feature.get("type") not in UNIPROT_FEATURE_TYPES:
            continue
        location = feature.get("location", {})
        domains.append({
            "name": feature.get("description") or feature.get("type"),
            "start": (location.get("start") or {}).get("value"),
            "end": (location.get("end") or {}).get("value"),
        })
    return {"protein_length": length, "domains": domains}


@lru_cache(maxsize=2048)
def fetch_uniprot(accession):
    return parse_uniprot(http("GET", f"{UNIPROT_URL}/{accession}.json"))


def build_feature_record(gene, uniprot_info=None):
    """Turn a MyGene.info gene record (plus optional UniProt info) into one tidy dict."""
    position = _as_list(gene.get("genomic_pos"))
    position = position[0] if position else {}
    n_transcripts, max_exons = _exon_summary(gene.get("exons"))
    accession = _first_swissprot(gene.get("uniprot"))
    uniprot_info = uniprot_info or {}
    go = gene.get("go") or {}
    return {
        "symbol": gene.get("symbol"),
        "name": gene.get("name"),
        "summary": gene.get("summary"),
        "structural": {
            "map_location": gene.get("map_location"),
            "chromosome": position.get("chr"),
            "strand": "+" if position.get("strand") == 1 else "-" if position.get("strand") == -1 else None,
            "transcripts": n_transcripts or None,
            "max_exons": max_exons,
            "uniprot": accession,
            "protein_length": uniprot_info.get("protein_length"),
            "domains": uniprot_info.get("domains", []),
            "alphafold_url": f"{ALPHAFOLD_ENTRY_URL}/{accession}" if accession else None,
        },
        "functional": {
            "gene_type": gene.get("type_of_gene"),
            "go_process": _go_terms(go, "BP"),
            "go_function": _go_terms(go, "MF"),
            "go_component": _go_terms(go, "CC"),
        },
    }


def fetch_features(valid):
    """Structural and functional features for validated genes (dict from validate_symbols)."""
    ids = [info["entrez"] for info in valid.values()]
    records = []
    for start in range(0, len(ids), BATCH_SIZE):
        genes = http("POST", f"{MYGENE_URL}/gene", data={
            "ids": ",".join(ids[start:start + BATCH_SIZE]), "fields": GENE_FIELDS, "species": "human"})
        for gene in genes:
            if gene.get("notfound"):
                continue
            accession = _first_swissprot(gene.get("uniprot"))
            info = None
            if accession:
                try:
                    info = fetch_uniprot(accession)
                except ApiError:
                    info = None
            records.append(build_feature_record(gene, info))
    return records


# ======================================================== MODULE 2: pathways

def parse_enrichr(rows, top=TOP_PATHWAYS):
    """Enrichr row: [rank, term, p-value, z-score, combined score, genes, adj. p-value, ...]."""
    pathways = []
    for row in rows:
        term = str(row[1])
        match = re.search(r"\s(R-HSA-\d+)$", term)
        pathways.append({
            "pathway": term[:match.start()] if match else term,
            "reactome_id": match.group(1) if match else None,
            "p_value": float(row[2]),
            "adj_p_value": float(row[6]),
            "genes": sorted(str(g).upper() for g in row[5]),
        })
    significant = [p for p in pathways if p["adj_p_value"] < 0.05]
    chosen = significant or pathways
    chosen.sort(key=lambda p: (p["adj_p_value"], -len(p["genes"])))
    for p in chosen:
        p["gene_count"] = len(p["genes"])
        p["significant"] = p["adj_p_value"] < 0.05
    return chosen[:top]


def enrich_pathways(symbols):
    """Pathway enrichment with Enrichr's Reactome library (hypergeometric test, FDR corrected)."""
    added = http("POST", f"{ENRICHR_URL}/addList",
                 files={"list": (None, "\n".join(symbols)), "description": (None, "GeneDisease Finder")})
    list_id = added.get("userListId")
    if list_id is None:
        raise ApiError("Enrichr did not accept the gene list")
    for library in ENRICHR_LIBRARIES:
        result = http("GET", f"{ENRICHR_URL}/enrich", params={"userListId": list_id, "backgroundType": library})
        rows = result.get(library) or []
        if rows:
            return parse_enrichr(rows)
    return []


def common_features(features, pathways):
    """GO terms and pathways shared by two or more of the input genes."""
    shared = defaultdict(set)
    for gene in features:
        for branch, label in (("go_process", "GO process"), ("go_function", "GO function"), ("go_component", "GO component")):
            for term in gene["functional"][branch]:
                shared[(label, term)].add(gene["symbol"])
    for p in pathways:
        shared[("Pathway", p["pathway"])].update(p["genes"])
    rows = [{"type": t, "feature": f, "genes": sorted(g), "count": len(g)} for (t, f), g in shared.items() if len(g) >= 2]
    rows.sort(key=lambda r: (-r["count"], r["type"], r["feature"]))
    return rows[:20]


def fetch_string_network(symbols, min_score=0.4):
    """Protein-protein interactions between the input genes from STRING."""
    edges = http("POST", STRING_URL, data={
        "identifiers": "\r".join(symbols[:MAX_NETWORK_GENES]), "species": 9606,
        "required_score": int(min_score * 1000), "caller_identity": "genedisease_finder"})
    return parse_string(edges)


def parse_string(edges):
    seen, out = set(), []
    for e in edges:
        a, b = str(e.get("preferredName_A", "")).upper(), str(e.get("preferredName_B", "")).upper()
        if not a or not b or a == b:
            continue
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        out.append({"a": key[0], "b": key[1], "score": float(e.get("score", 0))})
    degree = defaultdict(int)
    for e in out:
        degree[e["a"]] += 1
        degree[e["b"]] += 1
    hubs = sorted(degree.items(), key=lambda kv: (-kv[1], kv[0]))[:5]
    return {"edges": out, "hubs": hubs}


# ================================================== MODULE 3: diseases, drugs

OT_QUERY = """
query($id: String!) {
  target(ensemblId: $id) {
    approvedSymbol
    associatedDiseases(page: {index: 0, size: 50}) { rows { score disease { id name } } }
    tractability { label modality value }
  }
}"""


def parse_target(data):
    """Open Targets response -> ({disease: (id, score)}, [druggability labels])."""
    target = ((data or {}).get("data") or {}).get("target") or {}
    diseases = {}
    for row in (target.get("associatedDiseases") or {}).get("rows", []):
        disease = row.get("disease") or {}
        if disease.get("name"):
            diseases[disease["name"]] = (disease.get("id"), float(row.get("score", 0)))
    labels = sorted({t["label"] for t in target.get("tractability") or []
                     if t.get("value") and t.get("label") in DRUGGABLE_LABELS})
    return diseases, labels


def _fetch_one_target(item):
    symbol, ensembl = item
    if not ensembl:
        return symbol, {}, []
    data = http("POST", OPENTARGETS_URL, json={"query": OT_QUERY, "variables": {"id": ensembl}})
    diseases, labels = parse_target(data)
    return symbol, diseases, labels


def fetch_gene_diseases(valid):
    """Diseases and druggability for each validated gene, fetched in parallel."""
    items = [(info["symbol"], info["ensembl"]) for info in valid.values()][:MAX_DISEASE_GENES]
    gene_diseases, druggable, failures = {}, {}, 0
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = [pool.submit(_fetch_one_target, item) for item in items]
        for future in futures:
            try:
                symbol, diseases, labels = future.result()
            except ApiError:
                failures += 1
                continue
            gene_diseases[symbol] = diseases
            if labels:
                druggable[symbol] = labels
    if items and failures == len(items):
        raise ApiError("api.platform.opentargets.org: no requests succeeded")
    return gene_diseases, druggable, failures


def rank_diseases(gene_diseases, min_score=0.0, hide_generic=False):
    """Group by disease; rank by number of supporting genes, then mean score."""
    grouped = defaultdict(dict)
    for gene, diseases in gene_diseases.items():
        for name, (_id, score) in diseases.items():
            if hide_generic and name.strip().lower() in GENERIC_DISEASES:
                continue
            if score >= min_score:
                grouped[name][gene] = score
    rows = []
    for name, genes in grouped.items():
        scores = list(genes.values())
        rows.append({
            "disease": name, "n_genes": len(genes), "genes": sorted(genes),
            "mean_score": round(sum(scores) / len(scores), 3), "max_score": round(max(scores), 3),
            "scores": {g: round(s, 3) for g, s in genes.items()},
        })
    rows.sort(key=lambda r: (-r["n_genes"], -r["mean_score"], r["disease"]))
    return rows


def related_diseases(ranked, top=5):
    """Diseases that share the most genes with the top-ranked disease."""
    if not ranked:
        return []
    lead = set(ranked[0]["genes"])
    rows = [{"disease": r["disease"], "shared": len(lead & set(r["genes"]))} for r in ranked[1:]]
    rows = [r for r in rows if r["shared"] > 0]
    rows.sort(key=lambda r: (-r["shared"], r["disease"]))
    return rows[:top]


OT_EXPR_QUERY = """
query($id: String!) {
  target(ensemblId: $id) {
    baselineExpression(page: {index: 0, size: 200}) {
      rows { tissueBiosample { biosampleName } median unit datasourceId }
    }
  }
}"""


def parse_expression(data, top=5):
    """Open Targets baseline expression -> the tissues with the highest median TPM (GTEx rows when present)."""
    if (data or {}).get("errors") and not (data or {}).get("data"):
        raise ApiError("api.platform.opentargets.org: " + str(data["errors"])[:150])
    target = ((data or {}).get("data") or {}).get("target") or {}
    rows = (target.get("baselineExpression") or {}).get("rows") or []
    gtex = [r for r in rows if "gtex" in str(r.get("datasourceId", "")).lower()]
    best = {}
    for row in gtex or rows:
        label = (row.get("tissueBiosample") or {}).get("biosampleName")
        value = row.get("median")
        if label and value is not None and float(value) > 0:
            best[label] = max(best.get(label, 0.0), float(value))
    ranked = sorted(best.items(), key=lambda kv: -kv[1])
    return [{"tissue": t, "tpm": round(v, 1)} for t, v in ranked[:top]]


def _fetch_one_expression(item):
    symbol, ensembl = item
    if not ensembl:
        return symbol, []
    data = http("POST", OPENTARGETS_URL, json={"query": OT_EXPR_QUERY, "variables": {"id": ensembl}})
    return symbol, parse_expression(data)


def fetch_expression(valid):
    """Top expressing tissues per gene (Open Targets baseline expression), fetched in parallel."""
    items = [(info["symbol"], info["ensembl"]) for info in valid.values()]
    out, failures, first_error = {}, 0, ""
    with ThreadPoolExecutor(max_workers=8) as pool:
        for future in [pool.submit(_fetch_one_expression, i) for i in items]:
            try:
                symbol, rows = future.result()
            except ApiError as exc:
                failures += 1
                first_error = first_error or str(exc)
                continue
            out[symbol] = rows
    if items and failures == len(items):
        raise ApiError(first_error)
    return out


def _count(query):
    """Number of ClinVar-annotated variants matching a MyVariant.info query."""
    data = http("GET", MYVARIANT_URL, params={"q": query, "size": 1, "fields": "_id"})
    return int(data.get("total", 0))


def parse_variants(total, pathogenic):
    return {"total": int(total), "pathogenic": int(pathogenic)}


def _fetch_one_variants(symbol):
    gene = f"clinvar.gene.symbol:{symbol}"
    total = _count(gene)
    pathogenic = 0
    if total:
        try:
            pathogenic = _count(f"{gene} AND clinvar.rcv.clinical_significance:pathogenic")
        except ApiError:
            pathogenic = _count(f"{gene} AND clinvar.clinical_significance:pathogenic")
    return symbol, parse_variants(total, pathogenic)


def fetch_variants(valid):
    """ClinVar variant counts per gene through MyVariant.info."""
    symbols = [info["symbol"] for info in valid.values()]
    out, failures, first_error = {}, 0, ""
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in [pool.submit(_fetch_one_variants, sym) for sym in symbols]:
            try:
                symbol, counts = future.result()
            except ApiError as exc:
                failures += 1
                first_error = first_error or str(exc)
                continue
            out[symbol] = counts
    if symbols and failures == len(symbols):
        raise ApiError(first_error)
    return out


def key_findings(r):
    """A few plain sentences that summarize the results."""
    out = []
    if r["pathways"]:
        p = r["pathways"][0]
        out.append(f"Top pathway: {p['pathway']} (adj. p = {p['adj_p_value']:.2e}) with {', '.join(p['genes'])}.")
    if r["ranked"]:
        d = r["ranked"][0]
        out.append(f"Top disease: {d['disease']}, supported by {d['n_genes']} of {len(r['valid_symbols'])} genes (mean score {d['mean_score']:.2f}).")
    if r.get("string") and r["string"]["hubs"]:
        gene, count = r["string"]["hubs"][0]
        out.append(f"Most connected gene: {gene} ({count} interaction{'s' if count != 1 else ''}).")
    approved = sorted(g for g, labels in r["druggable"].items() if "Approved Drug" in labels)
    if approved:
        out.append("Approved-drug targets: " + ", ".join(approved) + ".")
    if r.get("gene_diseases") is not None:
        linked = {g for d in r["ranked"] for g in d["genes"]}
        missing = [g for g in r["valid_symbols"] if g in r["gene_diseases"] and g not in linked]
        if missing:
            out.append("No disease link at the current filter: " + ", ".join(missing) + ".")
    if r.get("compare"):
        c = r["compare"]
        extra = f" and {len(c['shared_pathways'])} pathways" if "shared_pathways" in c else ""
        out.append(f"Sets A and B share {len(c['shared'])} genes{extra}.")
    return out


def refresh_disease_views(result):
    """Re-rank diseases for the current filters, then rebuild everything that depends on them."""
    result["ranked"] = rank_diseases(result.get("gene_diseases", {}), result["min_score"], result["hide_generic"])
    result["related"] = related_diseases(result["ranked"])
    result["charts"] = build_chart_data(result["pathways"], result["ranked"], result["string"])
    result["findings"] = key_findings(result)
    result["stats"] = {
        "pathways": len(result["pathways"]), "diseases": len(result["ranked"]),
        "interactions": len((result.get("string") or {}).get("edges", [])),
        "approved": sum(1 for labels in result["druggable"].values() if "Approved Drug" in labels),
    }


# ========================================== MODULE 4: compare, charts, reports

def compare_gene_sets(symbols_a, symbols_b, pathways_a=None, pathways_b=None):
    set_a, set_b = set(symbols_a), set(symbols_b)
    out = {"shared": sorted(set_a & set_b), "only_a": sorted(set_a - set_b), "only_b": sorted(set_b - set_a)}
    if pathways_a is not None and pathways_b is not None:
        names_a = {p["pathway"] for p in pathways_a}
        names_b = {p["pathway"] for p in pathways_b}
        out["shared_pathways"] = sorted(names_a & names_b)
    return out


def circle_layout(nodes):
    count = max(len(nodes), 1)
    return {n: (math.cos(2 * math.pi * i / count), math.sin(2 * math.pi * i / count)) for i, n in enumerate(nodes)}


def build_chart_data(pathways, ranked, string_net):
    """Everything the page's Plotly charts need, as plain lists."""
    chart = {"pathways": [], "diseases": [], "heatmap": None, "ppi": None, "gene_disease": None}
    chart["pathways"] = [{"label": p["pathway"], "value": p["gene_count"]} for p in pathways[:10]]
    top = ranked[:10]
    chart["diseases"] = [{"label": r["disease"], "value": r["n_genes"]} for r in top]
    if top:
        genes = sorted({g for r in top for g in r["genes"]})
        chart["heatmap"] = {
            "genes": genes, "diseases": [r["disease"] for r in top],
            "z": [[r["scores"].get(g) for g in genes] for r in top],
        }
        top5 = top[:6]
        gene_nodes = sorted({g for r in top5 for g in r["genes"]})[:30]
        chart["gene_disease"] = {
            "genes": gene_nodes, "diseases": [r["disease"] for r in top5],
            "edges": [[g, r["disease"]] for r in top5 for g in r["genes"] if g in gene_nodes],
        }
    if string_net and string_net["edges"]:
        nodes = sorted({n for e in string_net["edges"] for n in (e["a"], e["b"])})
        chart["ppi"] = {"pos": circle_layout(nodes), "edges": [[e["a"], e["b"]] for e in string_net["edges"]]}
    return chart


def csv_response(rows, header, filename):
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(header)
    writer.writerows(rows)
    return Response(buffer.getvalue(), mimetype="text/csv",
                    headers={"Content-Disposition": f"attachment; filename={filename}"})


# ============================================================ run everything

def analyze(symbols_a, symbols_b=None, min_score=0.0, hide_generic=False):
    """Run all four modules. A failing database adds a warning instead of stopping the page."""
    warnings = []
    valid, invalid, notes = validate_symbols(symbols_a)          # Module 1 (needed by everything else)
    official = [info["symbol"] for info in valid.values()]
    result = {
        "input_count": len(symbols_a), "valid_symbols": official,
        "chromosomes": {c: n for c in CHROMOSOMES if (n := sum(1 for v in valid.values() if v.get("chr") == c))}, "invalid": invalid, "notes": notes,
        "features": [], "features_limited": len(valid) > MAX_FEATURE_GENES,
        "pathways": [], "common": [], "string": None,
        "ranked": [], "related": [], "druggable": {}, "min_score": min_score, "hide_generic": hide_generic,
        "compare": None, "warnings": warnings,
    }
    if not valid:
        return result

    def attempt(label, func):
        try:
            return func()
        except (ApiError, KeyError, IndexError, TypeError, ValueError) as exc:
            warnings.append(f"{label} unavailable: {exc}")
            return None

    shown = dict(list(valid.items())[:MAX_FEATURE_GENES])
    result["features"] = attempt("Gene features", lambda: fetch_features(shown)) or []
    expression = attempt("Tissue expression", lambda: fetch_expression(shown)) or {}
    variants = attempt("Variants", lambda: fetch_variants(shown)) or {}
    for record in result["features"]:
        record["expression"] = expression.get(record["symbol"], [])
        record["variants"] = variants.get(record["symbol"])
    result["pathways"] = attempt("Pathways", lambda: enrich_pathways(official)) or []
    result["common"] = common_features(result["features"], result["pathways"])
    result["string"] = attempt("Interaction network", lambda: fetch_string_network(official))

    fetched = attempt("Diseases", lambda: fetch_gene_diseases(valid))
    if fetched:
        gene_diseases, result["druggable"], failures = fetched
        result["gene_diseases"] = gene_diseases
        if failures:
            warnings.append(f"Diseases: {failures} gene lookups failed and were skipped.")

    if symbols_b:
        valid_b, invalid_b, _ = validate_symbols(symbols_b)
        official_b = [info["symbol"] for info in valid_b.values()]
        pathways_b = attempt("Pathways for set B", lambda: enrich_pathways(official_b)) if official_b else None
        result["compare"] = compare_gene_sets(official, official_b, result["pathways"], pathways_b)
        result["compare"]["invalid_b"] = [i["symbol"] for i in invalid_b]
        result["compare"]["count_a"], result["compare"]["count_b"] = len(official), len(official_b)

    refresh_disease_views(result)
    return result


# ===================================================================== pages

PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>GeneDisease Finder</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,500;12..96,700&family=Public+Sans:wght@400;500;600&display=swap" rel="stylesheet">
<style>
:root{--ink:#1a1147;--hero:#1a1147;--teal:#00a6c8;--coral:#ff3d81;--amber:#ffb000;--violet:#7a3cff;--lime:#14c38e;
--bg:#f6f3ff;--paper:#fff;--line:#e2dcf5;--muted:#625b82;--head:'Bricolage Grotesque','Segoe UI',Arial,sans-serif;--body:'Public Sans','Segoe UI',Arial,sans-serif}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:84px}
body{margin:0;font-family:var(--body);font-size:16px;line-height:1.55;background:var(--bg);color:var(--ink)}
a{color:var(--violet)}h1,h2,h3{font-family:var(--head);margin:0;letter-spacing:-.015em;line-height:1.15}
:focus-visible{outline:3px solid var(--amber);outline-offset:2px}
.top{position:sticky;top:0;z-index:20;background:linear-gradient(90deg,#1a1147,#3b1a8f);color:#fff;display:flex;align-items:center;justify-content:space-between;padding:10px 24px}
.brand{display:flex;align-items:center;gap:10px;color:#fff;text-decoration:none;font:700 20px var(--head)}
.brand svg{width:26px;height:26px}.top .quiet{color:#b9c9dc;font-size:14px}
.hero{background:radial-gradient(700px 300px at 85% 0,rgba(255,61,129,.45),transparent 70%),radial-gradient(600px 300px at 10% 100%,rgba(0,166,200,.4),transparent 70%),linear-gradient(135deg,#1a1147,#4a1fb0);color:#fff;padding:36px 24px 8px}
.wrap{max-width:1160px;margin:0 auto}
.hero h1{font-size:clamp(30px,4.4vw,46px);font-weight:700;max-width:16ch}
.hero .sub{color:#d9d0ff;margin:8px 0 0;max-width:60ch}
.ledger{display:flex;flex-wrap:wrap;gap:0;margin:22px 0 0;border-top:1px solid rgba(255,255,255,.25)}
.ledger div{padding:14px 28px 12px 0;margin-right:28px;border-right:1px solid rgba(255,255,255,.2)}.ledger div:last-child{border-right:0}
.ledger b{display:block;font:700 34px/1 var(--head);font-variant-numeric:tabular-nums}.ledger span{color:#d9d0ff;font-size:14px}
.ledger div:nth-child(1) b{color:#5ee7ff}.ledger div:nth-child(2) b{color:#ff7eb3}.ledger div:nth-child(3) b{color:#c3a6ff}.ledger div:nth-child(4) b{color:#ffd166}
.ideo{display:block;width:100%;max-height:230px;margin-top:6px}
.ideo .chr{transform-box:fill-box;transform-origin:50% 100%;animation:rise .7s cubic-bezier(.2,.8,.2,1) both}
@keyframes rise{from{transform:scaleY(.04);opacity:0}to{transform:scaleY(1);opacity:1}}
@media (prefers-reduced-motion:reduce){.ideo .chr{animation:none}html{scroll-behavior:auto}}
.cap{color:#bdb0f0;font-size:13px;margin:0 0 20px}
.layout{max-width:1160px;margin:0 auto;padding:28px 24px 60px;display:grid;grid-template-columns:230px minmax(0,1fr);gap:40px}
.rail{position:sticky;top:76px;align-self:start;font-size:15px}
.rail nav a{display:flex;justify-content:space-between;gap:8px;padding:7px 10px;margin:0 -10px;border-radius:8px;color:var(--ink);text-decoration:none}
.rail nav a:before{content:'';width:10px;height:10px;border-radius:50%;background:var(--dot,#7a3cff);align-self:center;flex:none}.rail nav a{align-items:center;justify-content:flex-start}.rail nav a small{margin-left:auto}.rail nav a:hover{background:#ece6ff}.rail nav a.on{background:var(--ink);color:#fff}.rail nav a small{opacity:.7;font-variant-numeric:tabular-nums}
.rail h3{font-size:15px;margin:22px 0 8px}.rail .btn{display:block;margin:0 0 6px;text-align:left}
.filter{background:var(--paper);border:1px solid var(--line);border-top:4px solid var(--amber);border-radius:12px;padding:12px 14px;margin-top:22px}
.filter h3{margin:0 0 6px}.filter label{display:block;font-size:14px;margin:6px 0 0}.filter input[type=range]{width:100%;accent-color:var(--coral)}
main{min-width:0;max-width:860px}
section{margin:0 0 54px}section>h2{font-size:30px;display:inline-block;background:linear-gradient(var(--sec,#7a3cff),var(--sec,#7a3cff)) 0 100%/100% 8px no-repeat;padding-bottom:2px}
#features{--sec:rgba(122,60,255,.35)}#pathways{--sec:rgba(0,166,200,.4)}#diseases{--sec:rgba(255,61,129,.35)}#compare{--sec:rgba(255,176,0,.45)}section>.lead{color:var(--muted);margin:6px 0 18px;max-width:62ch}
.findings{background:linear-gradient(135deg,#fff,#fbf7ff);border:1px solid var(--line);border-left:6px solid var(--coral);border-radius:14px;padding:16px 22px;margin-bottom:46px}
.findings h2{font-size:22px}.findings ul{list-style:none;margin:10px 0 0;padding:0}
.findings li{padding:9px 0 9px 22px;border-bottom:1px solid var(--line);position:relative;font-size:17px}
.findings li:before{content:"";position:absolute;left:0;top:19px;width:10px;height:10px;border-radius:50%;background:linear-gradient(135deg,var(--coral),var(--violet))}
.panel{background:var(--paper);border:1px solid var(--line);border-top:4px solid var(--violet);border-radius:14px;padding:18px 22px;margin-bottom:18px}
.panel h3{font-size:18px;margin-bottom:8px}
.notice{border-radius:10px;padding:10px 14px;margin:0 0 10px;font-size:15px}
.alert{background:#fdeceb;border:1px solid #f1b8b4;color:#8a1f17}.warn{background:#fff6df;border:1px solid #ecd49a;color:#6b4e00}.info{background:#e6eff4;border:1px solid #cfdde6}
.muted{color:var(--muted)}.error{color:#a32d2d}
.specs{background:var(--paper);border:1px solid var(--line);border-radius:14px;overflow:hidden}
.spec{border-bottom:1px solid var(--line)}.spec:last-child{border-bottom:0}
.spec>summary{list-style:none;cursor:pointer;display:flex;flex-wrap:wrap;align-items:center;gap:6px 10px;padding:12px 18px}
.spec>summary::-webkit-details-marker{display:none}
.spec>summary:before{content:"";width:8px;height:8px;border-right:2px solid var(--muted);border-bottom:2px solid var(--muted);transform:rotate(-45deg);transition:transform .15s}
.spec[open]>summary:before{transform:rotate(45deg)}.spec[open]>summary{background:#f3eeff}
.sym{font:700 18px var(--head);min-width:76px}.nm{color:var(--muted);flex:1 1 200px;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.chip{display:inline-block;background:#d8f6fc;color:#005f73;border-radius:6px;padding:1px 8px;font-size:13px;white-space:nowrap}
.chip.v{background:#ece3ff;color:#4b1fb0}.chip.d{background:#ffe9a8;color:#6b4300}.chip.b{background:#ffe0ec;color:#a3124f}.chip.bad{background:#fdeceb;color:#8a1f17}
.body{padding:6px 22px 20px}.two{display:grid;grid-template-columns:1fr 1fr;gap:26px}
h4{font:700 16px var(--head);margin:14px 0 6px}h4.s{color:var(--violet)}h4.f{color:var(--coral)}
table{width:100%;border-collapse:collapse;font-size:15px}
th{text-align:left;color:var(--muted);font-weight:400;vertical-align:top;padding:5px 10px 5px 0;border-top:1px solid var(--line)}
td{padding:5px 8px 5px 0;border-top:1px solid var(--line);vertical-align:top}table.grid th{font-weight:600;color:var(--ink)}
.meter{display:inline-block;vertical-align:middle;height:8px;border-radius:4px;background:linear-gradient(90deg,#00a6c8,#7a3cff);margin-right:8px}.meter.c{background:linear-gradient(90deg,#ff3d81,#ffb000)}
.chart{width:100%;height:380px}.vsplit{display:grid;gap:18px}
.btn,button{padding:9px 16px;border:2px solid var(--violet);border-radius:9px;background:#fff;color:var(--ink);font:600 15px var(--body);cursor:pointer;text-decoration:none;display:inline-block}
button.primary{background:linear-gradient(135deg,#7a3cff,#ff3d81);border-color:transparent;color:#fff}.btn:hover,button:hover{background:#ece6ff}button.primary:hover{filter:brightness(1.1);background:linear-gradient(135deg,#7a3cff,#ff3d81)}
textarea,input[type=file],input[type=number]{width:100%;padding:10px;border:1px solid var(--line);border-radius:9px;background:#fff;font:inherit}
label.f{display:block;margin:16px 0 6px;font-weight:600}
.sources{margin-top:8px;color:var(--muted);font-size:14px}.sources .chip{background:#ece6ff;color:var(--ink)}
.intro{background:radial-gradient(700px 360px at 90% 0,rgba(255,61,129,.45),transparent 70%),radial-gradient(600px 320px at 0 100%,rgba(0,166,200,.4),transparent 70%),linear-gradient(135deg,#1a1147,#4a1fb0);color:#fff;padding:54px 24px 0}
.intro .grid{max-width:1160px;margin:0 auto;display:grid;grid-template-columns:minmax(0,1fr) minmax(0,480px);gap:48px;align-items:start}
.intro h1{font-size:clamp(34px,5vw,54px);font-weight:700;max-width:15ch}.intro p{color:#d9d0ff;max-width:46ch;margin:14px 0 0;font-size:18px}
.intro ul{list-style:none;padding:0;margin:22px 0 0;color:#ece6ff}.intro li{padding:6px 0 6px 20px;position:relative}.intro li:before{content:'';position:absolute;left:0;top:15px;width:10px;height:10px;border-radius:50%;background:var(--c)}.intro li:nth-child(1){--c:#b794ff}.intro li:nth-child(2){--c:#5ee7ff}.intro li:nth-child(3){--c:#ff7eb3}.intro li:nth-child(4){--c:#ffd166}.intro li b{color:#fff}
.formcard{background:var(--paper);color:var(--ink);border-radius:16px;border-top:6px solid var(--amber);box-shadow:0 20px 50px rgba(10,0,60,.35);padding:22px 24px 24px;margin-bottom:34px}
.formcard h2{font-size:22px}.formcard .buttons{margin-top:18px;display:flex;gap:10px;flex-wrap:wrap}
@media(max-width:900px){.layout{grid-template-columns:1fr;gap:12px}.rail{position:static}.rail nav{display:flex;flex-wrap:wrap;gap:4px}.rail nav a{margin:0}
.intro .grid{grid-template-columns:1fr;gap:10px}.two{grid-template-columns:1fr}.ledger div{padding-right:18px;margin-right:18px}}
@page{margin:12mm}
@media print{*{-webkit-print-color-adjust:exact;print-color-adjust:exact}.top,.rail,.noprint{display:none}.layout{display:block;padding:12px 0}
html{scroll-behavior:auto}.spec>summary:before{display:none}details{display:block}details>*{display:block}.spec>summary{display:flex}
.chart,tr,.chip,svg{break-inside:avoid;page-break-inside:avoid}h2,h3{break-after:avoid;page-break-after:avoid}.ideo .chr{animation:none}}
</style></head><body>
<header class="top"><a class="brand" href="{{ url_for('index') }}"><svg viewBox="0 0 26 26" aria-hidden="true"><path d="M5 2c0 8 16 9 16 18M21 2c0 8-16 9-16 18" fill="none" stroke="#5ee7ff" stroke-width="2.4" stroke-linecap="round"/><path d="M7 6h12M6 11h14M7 16h12" stroke="#ff7eb3" stroke-width="1.8" stroke-linecap="round"/></svg>GeneDisease Finder</a>
<span class="quiet">{% if r is none %}Gene list in, diseases out{% else %}Results{% endif %}</span></header>

{% if r is none %}
<section class="intro"><div class="grid">
<div><h1>Find the diseases behind a gene list</h1>
<p>Paste or upload your genes. You get their structure and function, the pathways they share, and the diseases they point to.</p>
<ul><li><b>Features.</b> Location, domains, variants, tissues</li><li><b>Pathways.</b> Where your genes work together</li><li><b>Diseases.</b> Ranked by evidence, with drug targets</li><li><b>Compare.</b> Set A against set B</li></ul></div>
<div class="formcard">
{% with messages = get_flashed_messages() %}{% for m in messages %}<p class="notice alert">{{ m }}</p>{% endfor %}{% endwith %}
<h2>Your genes</h2>
<form method="post" action="{{ url_for('run') }}" enctype="multipart/form-data">
<label class="f" for="genes">Gene set A (commas, spaces or new lines)</label>
<textarea id="genes" name="genes" rows="5" placeholder="TP53, BRCA1, EGFR, KRAS, PTEN"></textarea>
<label class="f" for="file">Or upload a file for set A (.csv, .tsv, .txt, .xlsx, up to 5 MB)</label>
<input id="file" type="file" name="file" accept=".csv,.tsv,.txt,.xlsx,.xls">
<label class="f" for="genes_b">Gene set B (optional, to compare with set A)</label>
<textarea id="genes_b" name="genes_b" rows="3" placeholder="IL6, TNF, INS"></textarea>
<label class="f" for="min_score">Minimum disease evidence score (0 to 1)</label>
<input id="min_score" type="number" name="min_score" min="0" max="1" step="0.05" value="0.2" style="max-width:140px">
<label class="f" style="font-weight:400"><input type="checkbox" name="hide_generic" value="1"> Hide very generic diseases (cancer, neoplasm, carcinoma...)</label>
<div class="buttons"><button class="primary" type="submit">Analyze genes</button>
<button type="submit" name="example" value="1">Load example genes</button></div></form></div></div>
<div class="wrap"><svg class="ideo" id="ideo" aria-hidden="true"></svg></div></section>
{% else %}
<section class="hero"><div class="wrap">
<h1>{{ r.valid_symbols|length }} of {{ r.input_count }} genes recognised</h1>
<p class="sub">Each bar below is a human chromosome, drawn to its relative length. Brighter bars hold more of your genes.</p>
<div class="ledger">
<div><b>{{ r.stats.pathways }}</b><span>pathways found</span></div>
<div><b>{{ r.stats.diseases }}</b><span>diseases linked</span></div>
<div><b>{{ r.stats.interactions }}</b><span>protein interactions</span></div>
<div><b>{{ r.stats.approved }}</b><span>approved-drug targets</span></div></div>
<svg class="ideo" id="ideo" role="img" aria-label="Number of your genes on each chromosome"></svg>
<p class="cap">Chromosome lengths are approximate. Genes without a known chromosome are not drawn.</p></div></section>

<div class="layout">
<aside class="rail noprint"><nav aria-label="Sections">
<a href="#features" style="--dot:#7a3cff">Gene features <small>{{ r.features|length }}</small></a>
<a href="#pathways" style="--dot:#00a6c8">Pathways <small>{{ r.pathways|length }}</small></a>
<a href="#diseases" style="--dot:#ff3d81">Diseases <small>{{ r.ranked|length }}</small></a>
{% if r.compare %}<a href="#compare" style="--dot:#ffb000">Compare sets</a>{% endif %}</nav>
{% if r.get('gene_diseases') is not none %}<form method="get" class="filter" action="{{ url_for('show', rid=rid) }}">
<h3>Disease filter</h3>
<label for="ms">Minimum evidence score <strong id="msv">{{ "%.2f"|format(r.min_score) }}</strong></label>
<input id="ms" type="range" name="min_score" min="0" max="0.9" step="0.05" value="{{ r.min_score }}" oninput="document.getElementById('msv').textContent=parseFloat(this.value).toFixed(2)" onchange="this.form.submit()">
<label><input type="checkbox" name="hide_generic" value="1" {% if r.hide_generic %}checked{% endif %} onchange="this.form.submit()"> Hide very generic diseases (cancer, neoplasm, carcinoma...)</label></form>{% endif %}
<h3>Export</h3>
<a class="btn" href="{{ url_for('download', rid=rid, kind='features') }}">Features CSV</a>
<a class="btn" href="{{ url_for('download', rid=rid, kind='pathways') }}">Pathways CSV</a>
<a class="btn" href="{{ url_for('download', rid=rid, kind='diseases') }}">Diseases CSV</a>
<button onclick="window.print()" style="width:100%;text-align:left">Print or save as PDF</button>
<a class="btn" style="margin-top:14px" href="{{ url_for('index') }}">New analysis</a></aside>

<main>
{% with messages = get_flashed_messages() %}{% for m in messages %}<p class="notice alert">{{ m }}</p>{% endfor %}{% endwith %}
{% for n in r.notes %}<p class="notice info">{{ n }}</p>{% endfor %}
{% if r.invalid %}<p class="notice alert">Not recognized: {% for i in r.invalid %}{{ i.symbol }}{% if i.suggestion %} (did you mean {{ i.suggestion }}?){% endif %}{% if not loop.last %}, {% endif %}{% endfor %}</p>{% endif %}
{% for w in r.warnings %}<p class="notice warn">{{ w }}</p>{% endfor %}
{% if r.findings %}<div class="findings"><h2>Key findings</h2><ul>{% for f in r.findings %}<li>{{ f }}</li>{% endfor %}</ul></div>{% endif %}

<section id="features"><h2>Gene features</h2>
<p class="lead">Structure and function of each gene. Open a row for details.{% if r.features_limited %} Shown for the first {{ max_features }} valid genes.{% endif %}</p>
<div class="specs">{% for g in r.features %}<details class="spec" {% if loop.first %}open{% endif %}>
<summary><span class="sym">{{ g.symbol }}</span><span class="nm">{{ g.name }}</span>
{% if g.structural.map_location %}<span class="chip">{{ g.structural.map_location }}</span>{% endif %}
{% if g.structural.protein_length %}<span class="chip">{{ g.structural.protein_length }} aa</span>{% endif %}
{% if g.variants and g.variants.pathogenic %}<span class="chip v">{{ g.variants.pathogenic }} pathogenic variants</span>{% endif %}
{% if 'Approved Drug' in r.druggable.get(g.symbol, []) %}<span class="chip d">Approved drug</span>{% endif %}</summary>
<div class="body">{% if g.summary %}<p>{{ g.summary }}</p>{% endif %}
<div class="two"><div><h4 class="s">Structural</h4><table>
<tr><th>Location</th><td>{{ g.structural.map_location or "-" }}</td></tr>
<tr><th>Chromosome</th><td>{{ g.structural.chromosome or "-" }}{% if g.structural.strand %} ({{ g.structural.strand }} strand){% endif %}</td></tr>
<tr><th>Transcripts</th><td>{{ g.structural.transcripts or "-" }}</td></tr>
<tr><th>Max exons</th><td>{{ g.structural.max_exons or "-" }}</td></tr>
<tr><th>Protein length</th><td>{% if g.structural.protein_length %}{{ g.structural.protein_length }} amino acids{% else %}-{% endif %}</td></tr>
<tr><th>Domains</th><td>{% for d in g.structural.domains %}{{ d.name }} ({{ d.start }}-{{ d.end }}){% if not loop.last %}<br>{% endif %}{% else %}-{% endfor %}</td></tr>
<tr><th>ClinVar variants</th><td>{% if g.variants and g.variants.total %}{{ g.variants.total }} in total, {{ g.variants.pathogenic }} pathogenic or likely pathogenic{% else %}-{% endif %}</td></tr>
<tr><th>3D model</th><td>{% if g.structural.alphafold_url %}<a href="{{ g.structural.alphafold_url }}" target="_blank" rel="noopener">AlphaFold DB</a>{% else %}-{% endif %}</td></tr></table></div>
<div><h4 class="f">Functional</h4><table>
<tr><th>Gene type</th><td>{{ g.functional.gene_type or "-" }}</td></tr>
<tr><th>GO process</th><td>{{ g.functional.go_process|join(", ") or "-" }}</td></tr>
<tr><th>GO function</th><td>{{ g.functional.go_function|join(", ") or "-" }}</td></tr>
<tr><th>GO component</th><td>{{ g.functional.go_component|join(", ") or "-" }}</td></tr>
<tr><th>Top tissues (RNA, TPM)</th><td>{% for t in g.expression %}{{ t.tissue }} ({{ t.tpm }}){% if not loop.last %}, {% endif %}{% else %}-{% endfor %}</td></tr>
<tr><th>Druggability</th><td>{% if r.druggable.get(g.symbol) %}{% for l in r.druggable[g.symbol] %}<span class="chip">{{ l }}</span> {% endfor %}{% else %}-{% endif %}</td></tr></table></div></div></div></details>{% endfor %}</div></section>

<section id="pathways"><h2>Pathways and common features</h2>
<p class="lead">Pathways where your genes appear more often than chance. The bar length is how many of your genes are in the pathway; rows are ordered by significance.</p>
{% if r.pathways %}{% set maxc = r.pathways|map(attribute='gene_count')|max %}
<div class="panel"><div id="c_path" class="chart"></div>
<table class="grid"><tr><th>Pathway</th><th>Your genes</th><th>Adj. p-value</th></tr>
{% for p in r.pathways %}<tr><td>{{ p.pathway }}</td><td><span class="meter" style="width:{{ (70 * p.gene_count / maxc)|round|int }}px"></span>{{ p.gene_count }}<br><span class="muted">{{ p.genes[:12]|join(", ") }}{% if p.genes|length > 12 %} and {{ p.genes|length - 12 }} more{% endif %}</span></td><td>{{ "%.2e"|format(p.adj_p_value) }}{% if not p.significant %} <span class="muted">(not significant)</span>{% endif %}</td></tr>{% endfor %}</table></div>
{% else %}<p class="muted">No pathways to show.</p>{% endif %}
{% if r.common %}<div class="panel"><h3>Common features (shared by 2 or more genes)</h3><table class="grid"><tr><th>Type</th><th>Feature</th><th>Genes</th></tr>
{% for c in r.common %}<tr><td>{{ c.type }}</td><td>{{ c.feature }}</td><td>{{ c.count }}: {{ c.genes[:12]|join(", ") }}{% if c.genes|length > 12 %} and {{ c.genes|length - 12 }} more{% endif %}</td></tr>{% endfor %}</table></div>{% endif %}
{% if r.string and r.string.edges %}<div class="panel"><h3>Interaction network (STRING)</h3><div id="c_ppi" class="chart"></div>
<p class="muted">A line joins two proteins that STRING reports as interacting. Position on the circle carries no meaning. Most connected genes: {% for g, n in r.string.hubs %}{{ g }} ({{ n }}){% if not loop.last %}, {% endif %}{% endfor %}.</p></div>{% endif %}</section>

<section id="diseases"><h2>Disease associations</h2>
<p class="lead">Diseases ranked by how many of your genes are linked to them, then by mean evidence score. Use the filter on the left to change the minimum score.</p>
{% if r.ranked %}<div class="panel"><div id="c_dis" class="chart"></div>
<table class="grid"><tr><th>#</th><th>Disease</th><th>Genes</th><th>Mean score</th></tr>
{% for d in r.ranked[:15] %}<tr><td>{{ loop.index }}</td><td>{{ d.disease }}</td><td>{{ d.n_genes }}: {{ d.genes[:12]|join(", ") }}{% if d.genes|length > 12 %} and {{ d.genes|length - 12 }} more{% endif %}</td><td><span class="meter c" style="width:{{ (50 * d.mean_score)|round|int }}px"></span>{{ "%.2f"|format(d.mean_score) }}</td></tr>{% endfor %}</table></div>
<div class="panel"><h3>Gene-disease heatmap (association score)</h3><div id="c_heat" class="chart" style="height:460px"></div></div>
<div class="panel"><h3>Gene-disease network</h3><div id="c_gd" class="chart" style="height:460px"></div></div>
{% if r.related %}<div class="panel"><h3>Related diseases (genes shared with {{ r.ranked[0].disease }})</h3>
{% for x in r.related %}<span class="chip b">{{ x.disease }} ({{ x.shared }})</span> {% endfor %}</div>{% endif %}
{% if r.druggable %}<div class="panel"><h3>Drug targets</h3><table class="grid"><tr><th>Gene</th><th>Druggability evidence</th></tr>
{% for g, labels in r.druggable.items() %}<tr><td>{{ g }}</td><td>{{ labels|join(", ") }}</td></tr>{% endfor %}</table></div>{% endif %}
{% else %}<p class="muted">No disease associations to show at this filter. Lower the minimum score on the left.</p>{% endif %}</section>

{% if r.compare %}<section id="compare"><h2>Compare gene sets</h2>
<p class="lead">Set A against set B: shared genes, genes found in only one set, and pathways both sets are enriched in.</p><div class="panel">
<p>Set A: <strong>{{ r.compare.count_a }}</strong> genes. Set B: <strong>{{ r.compare.count_b }}</strong> genes.
{% if r.compare.invalid_b %}<span class="error">Not recognized in B: {{ r.compare.invalid_b|join(", ") }}</span>{% endif %}</p>
<svg viewBox="0 0 600 240" width="100%" style="max-width:600px;display:block;margin:0 auto" role="img" aria-label="Genes only in set A, shared, and only in set B">
<circle cx="240" cy="110" r="95" fill="rgba(0,166,200,0.3)" stroke="#00a6c8"/><circle cx="360" cy="110" r="95" fill="rgba(255,61,129,0.3)" stroke="#ff3d81"/>
<g text-anchor="middle" fill="#0b2545" font-family="Public Sans, Arial">
<text x="200" y="108" font-size="28" font-weight="bold">{{ r.compare.only_a|length }}</text><text x="200" y="130" font-size="14">only A</text>
<text x="300" y="108" font-size="28" font-weight="bold">{{ r.compare.shared|length }}</text><text x="300" y="130" font-size="14">shared</text>
<text x="400" y="108" font-size="28" font-weight="bold">{{ r.compare.only_b|length }}</text><text x="400" y="130" font-size="14">only B</text>
<text x="200" y="228" font-size="15" font-weight="bold">Set A</text><text x="400" y="228" font-size="15" font-weight="bold">Set B</text></g></svg>
<table class="grid"><tr><th>Shared genes</th><td>{{ r.compare.shared|join(", ") or "none" }}</td></tr>
<tr><th>Only in A</th><td>{{ r.compare.only_a|join(", ") or "none" }}</td></tr>
<tr><th>Only in B</th><td>{{ r.compare.only_b|join(", ") or "none" }}</td></tr>
{% if r.compare.shared_pathways is defined %}<tr><th>Shared pathways</th><td>{{ r.compare.shared_pathways|join(", ") or "none" }}</td></tr>{% endif %}</table></div></section>{% endif %}

<p class="sources">Data from <span class="chip">MyGene.info</span> <span class="chip">UniProt</span> <span class="chip">Enrichr</span> <span class="chip">STRING</span> <span class="chip">Open Targets</span> <span class="chip">MyVariant.info</span>. These are associations from public databases, not diagnoses.</p>
</main></div>
{% endif %}

<script>
const CHROM = {{ chrom|tojson }};
const LEN = {"1":248,"2":242,"3":198,"4":190,"5":181,"6":171,"7":159,"8":145,"9":138,"10":134,"11":135,"12":133,"13":114,"14":107,"15":102,"16":90,"17":83,"18":80,"19":59,"20":64,"21":47,"22":51,"X":156,"Y":57};
const ACRO = ["13","14","15","21","22"];
function mix(a,b,t){const p=h=>[1,3,5].map(i=>parseInt(h.slice(i,i+2),16));const x=p(a),y=p(b);return 'rgb('+x.map((v,i)=>Math.round(v+(y[i]-v)*t)).join(',')+')';}
(function(){
  const el=document.getElementById('ideo'); if(!el) return;
  const keys=Object.keys(LEN), W=1000, H=210, bw=30, gap=(W-keys.length*bw)/(keys.length-1), top=170, sc=135/248;
  const maxc=Math.max(1,...keys.map(k=>CHROM[k]||0));
  let svg='';
  keys.forEach((k,i)=>{
    const x=i*(bw+gap), h=LEN[k]*sc, y=top-h, n=CHROM[k]||0, c=n?mix('#5b45b8','#ff7eb3',Math.sqrt(n/maxc)):'#5b45b8';
    const cy=y+h*(ACRO.includes(k)?0.12:0.4);
    svg+='<g class="chr" style="animation-delay:'+(i*32)+'ms"><title>Chromosome '+k+': '+n+(n===1?' gene':' genes')+'</title>'
      +'<rect x="'+x+'" y="'+y+'" width="'+bw+'" height="'+h+'" rx="15" fill="'+c+'"/>'
      +'<rect x="'+(x-1)+'" y="'+(cy-3)+'" width="'+(bw+2)+'" height="6" style="fill:var(--hero)"/>'
      +(n?'<text x="'+(x+bw/2)+'" y="'+(y-7)+'" text-anchor="middle" fill="#fff" font-size="13" font-weight="700">'+n+'</text>':'')+'</g>'
      +'<text x="'+(x+bw/2)+'" y="192" text-anchor="middle" fill="#8fa8c6" font-size="12">'+k+'</text>';
  });
  el.setAttribute('viewBox','0 0 '+W+' '+H); el.innerHTML=svg;
})();
(function(){
  const links=[...document.querySelectorAll('.rail nav a')]; if(!links.length||!('IntersectionObserver' in window)) return;
  const map=new Map(links.map(a=>[a.getAttribute('href').slice(1),a]));
  const io=new IntersectionObserver(es=>{es.forEach(e=>{if(e.isIntersecting){links.forEach(l=>l.classList.remove('on'));const a=map.get(e.target.id);if(a)a.classList.add('on');}});},{rootMargin:'-15% 0px -70% 0px'});
  map.forEach((a,id)=>{const s=document.getElementById(id);if(s)io.observe(s);});
})();
</script>
{% if r is not none %}
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script>
const D = {{ charts|tojson }};
const cfg = {displayModeBar:false, responsive:true};
const lay = (extra) => Object.assign({margin:{l:10,r:20,t:30,b:40}, paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)', font:{family:'Public Sans, Segoe UI, Arial', color:'#1a1147'}}, extra||{});
const short = (s) => s.length>56 ? s.slice(0,54)+'…' : s;
function bar(id, items, color, title, xtitle){
  const el=document.getElementById(id); if(!el||!items||!items.length) return;
  const top=Math.max(...items.map(i=>i.value));
  Plotly.newPlot(el,[{type:'bar',orientation:'h',x:items.map(i=>i.value).reverse(),y:items.map(i=>short(i.label)).reverse(),hovertext:items.map(i=>i.label).reverse(),hoverinfo:'text+x',marker:{color:color},text:items.map(i=>i.value).reverse(),textposition:'outside',cliponaxis:false}],
    lay({title:{text:title,font:{size:15}},margin:{l:320,r:40,t:40,b:40},xaxis:Object.assign({title:xtitle,rangemode:'tozero',tickformat:'d'},top<=10?{dtick:1}:{})}),cfg);
}
bar('c_path', D.pathways, '#00a6c8', 'Input genes in each enriched pathway', 'Number of input genes in the pathway');
bar('c_dis', D.diseases, '#ff3d81', 'Genes supporting each disease', 'Number of input genes linked to the disease');
if (D.heatmap && document.getElementById('c_heat')) {
  Plotly.newPlot('c_heat',[{type:'heatmap',x:D.heatmap.genes,y:D.heatmap.diseases,z:D.heatmap.z,colorscale:[[0,'#f3ecff'],[0.5,'#ff7eb3'],[1,'#7a3cff']],zmin:0,zmax:1,hoverongaps:false,colorbar:{title:'Score'}}],
    lay({margin:{l:260,r:20,t:20,b:80},xaxis:{tickangle:-45}}),cfg);
}
if (D.gene_disease && document.getElementById('c_gd')) {
  const g=D.gene_disease, gy={}, dy={};
  g.genes.forEach((n,i)=>gy[n]=-i/Math.max(g.genes.length-1,1)); g.diseases.forEach((n,i)=>dy[n]=-i/Math.max(g.diseases.length-1,1));
  const lx=[], ly=[]; g.edges.forEach(e=>{lx.push(0,1,null); ly.push(gy[e[0]],dy[e[1]],null);});
  Plotly.newPlot('c_gd',[
    {type:'scatter',mode:'lines',x:lx,y:ly,line:{color:'#b8c4cb',width:1},hoverinfo:'skip'},
    {type:'scatter',mode:'markers+text',x:g.genes.map(()=>0),y:g.genes.map(n=>gy[n]),text:g.genes,textposition:'middle left',cliponaxis:false,marker:{size:12,color:'#7a3cff'},hoverinfo:'text'},
    {type:'scatter',mode:'markers+text',x:g.diseases.map(()=>1),y:g.diseases.map(n=>dy[n]),text:g.diseases,textposition:'middle right',cliponaxis:false,marker:{size:14,color:'#ff3d81'},hoverinfo:'text'}],
    lay({showlegend:false,margin:{l:90,r:260,t:10,b:10},xaxis:{visible:false,range:[-0.2,1.2]},yaxis:{visible:false,range:[-1.1,0.1]}}),cfg);
}
if (D.ppi && document.getElementById('c_ppi')) {
  const p=D.ppi, ex=[], ey=[]; p.edges.forEach(e=>{ex.push(p.pos[e[0]][0],p.pos[e[1]][0],null); ey.push(p.pos[e[0]][1],p.pos[e[1]][1],null);});
  const nodes=Object.keys(p.pos);
  Plotly.newPlot('c_ppi',[
    {type:'scatter',mode:'lines',x:ex,y:ey,line:{color:'#b8c4cb',width:1},hoverinfo:'skip'},
    {type:'scatter',mode:'markers+text',x:nodes.map(n=>p.pos[n][0]),y:nodes.map(n=>p.pos[n][1]),text:nodes,textposition:'top center',cliponaxis:false,marker:{size:12,color:'#7a3cff'},hoverinfo:'text'}],
    lay({showlegend:false,margin:{l:20,r:20,t:20,b:20},xaxis:{visible:false},yaxis:{visible:false,scaleanchor:'x'}}),cfg);
}
</script>
{% endif %}
</body></html>"""


def render_page(result=None, rid=None):
    return render_template_string(
        PAGE, r=result, rid=rid, max_features=MAX_FEATURE_GENES,
        charts=(result or {}).get("charts") or {}, chrom=(result or {}).get("chromosomes") or {})


def store(result):
    rid = uuid.uuid4().hex[:10]
    RESULTS[rid] = result
    while len(RESULTS) > SAVED_RESULTS:
        RESULTS.popitem(last=False)
    return rid


@app.get("/")
def index():
    return render_page()


@app.post("/run")
def run():
    try:
        if request.form.get("example"):
            symbols_a, symbols_b = parse_genes_text(SAMPLE_GENES), parse_genes_text(SAMPLE_GENES_B)
        else:
            symbols_a = clean_symbols(parse_genes_text(request.form.get("genes", "")) + parse_genes_file(request.files.get("file")))
            symbols_b = parse_genes_text(request.form.get("genes_b", ""))
    except ValueError as exc:
        flash(str(exc))
        return redirect(url_for("index"))
    except Exception:
        flash("Could not read the file. Check that it is a valid CSV, TXT or Excel file.")
        return redirect(url_for("index"))
    if not symbols_a:
        flash("Enter at least one gene symbol or upload a file.")
        return redirect(url_for("index"))
    if len(symbols_a) > MAX_GENES or len(symbols_b) > MAX_GENES:
        flash(f"Too many genes. The limit is {MAX_GENES} per set.")
        return redirect(url_for("index"))
    try:
        min_score = min(max(float(request.form.get("min_score", 0.2)), 0.0), 1.0)
    except ValueError:
        min_score = 0.2
    try:
        result = analyze(symbols_a, symbols_b, min_score, bool(request.form.get("hide_generic")))
    except ApiError as exc:
        flash(f"Could not reach a database: {exc}")
        return redirect(url_for("index"))
    result["symbols_a"], result["symbols_b"] = symbols_a, symbols_b
    rid = store(result)
    return redirect(url_for("show", rid=rid))


@app.get("/results/<rid>")
def show(rid):
    result = RESULTS.get(rid)
    if result is None:
        flash("That result has expired. Run the analysis again.")
        return redirect(url_for("index"))
    if "min_score" in request.args:           # evidence score slider: re-rank without calling the APIs again
        try:
            score = min(max(float(request.args["min_score"]), 0.0), 1.0)
        except ValueError:
            score = result["min_score"]
        result["min_score"] = score
        result["hide_generic"] = "hide_generic" in request.args
        refresh_disease_views(result)
    return render_page(result, rid)


@app.get("/download/<rid>/<kind>.csv")
def download(rid, kind):
    result = RESULTS.get(rid)
    if result is None:
        abort(404)
    if kind == "features":
        rows = []
        for g in result["features"]:
            s, f = g["structural"], g["functional"]
            rows.append([g["symbol"], g["name"], s["map_location"], s["chromosome"], s["transcripts"], s["max_exons"],
                         s["protein_length"], "; ".join(f"{d['name']} ({d['start']}-{d['end']})" for d in s["domains"]),
                         f["gene_type"], "; ".join(f["go_process"]), "; ".join(f["go_function"]), "; ".join(f["go_component"]),
                         "; ".join(f"{t['tissue']} ({t['tpm']} TPM)" for t in g.get("expression", [])),
                         (g.get("variants") or {}).get("total", ""), (g.get("variants") or {}).get("pathogenic", "")])
        header = ["gene", "name", "location", "chromosome", "transcripts", "max_exons", "protein_length", "domains",
                  "gene_type", "go_process", "go_function", "go_component", "top_tissues", "clinvar_variants", "clinvar_pathogenic"]
    elif kind == "pathways":
        rows = [[p["pathway"], p["reactome_id"], p["gene_count"], ";".join(p["genes"]), p["p_value"], p["adj_p_value"]]
                for p in result["pathways"]]
        header = ["pathway", "reactome_id", "gene_count", "genes", "p_value", "adj_p_value"]
    elif kind == "diseases":
        rows = [[d["disease"], d["n_genes"], ";".join(d["genes"]), d["mean_score"], d["max_score"]] for d in result["ranked"]]
        header = ["disease", "n_genes", "genes", "mean_score", "max_score"]
    else:
        abort(404)
    return csv_response(rows, header, f"genedisease_{kind}.csv")


@app.errorhandler(413)
def too_large(_error):
    flash("File is too large. The limit is 5 MB.")
    return redirect(url_for("index"))


if __name__ == "__main__":
    # Local use: python app.py. Set HOST=0.0.0.0 to open it from a phone on the same Wi-Fi.
    app.run(host=os.environ.get("HOST", "127.0.0.1"), port=int(os.environ.get("PORT", "5000")), debug=False)
