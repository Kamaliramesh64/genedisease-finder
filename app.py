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
        hits = http("POST", f"{MYGENE_URL}/query", data={
            "q": ",".join(batch), "scopes": "symbol,alias", "species": "human",
            "fields": "symbol,entrezgene,ensembl.gene", "size": 1,
        })
        found = set()
        for hit in hits:
            query = str(hit.get("query", "")).upper()
            if hit.get("notfound") or "entrezgene" not in hit or query in found:
                continue
            found.add(query)
            official = str(hit.get("symbol", query)).upper()
            valid[query] = {"symbol": official, "entrez": str(hit["entrezgene"]), "ensembl": _ensembl_id(hit)}
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
        "input_count": len(symbols_a), "valid_symbols": official, "invalid": invalid, "notes": notes,
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
<style>
:root{--navy:#0b2545;--teal:#0f9488;--coral:#d4552f;--bg:#f4f8fa;--line:#d5dee3;--muted:#5b6770}
*{box-sizing:border-box}body{margin:0;font-family:Calibri,Arial,sans-serif;background:var(--bg);color:var(--navy);line-height:1.5}
header{background:var(--navy);padding:14px 24px}header a{color:#fff;font-weight:bold;font-size:20px;text-decoration:none}
main{max-width:1000px;margin:24px auto;padding:0 16px}h1{font-size:26px}h2{font-size:21px;margin-top:32px}
.muted{color:var(--muted)}.card{background:#fff;border:1px solid var(--line);border-radius:12px;padding:16px 20px;margin-bottom:14px}
.alert{background:#fdeceb;border:1px solid #f1b8b4;border-radius:8px;padding:10px 14px;color:#8a1f17}
.warn{background:#fff6df;border:1px solid #ecd49a;border-radius:8px;padding:8px 14px;color:#6b4e00;margin:6px 0}
.error{color:#a32d2d}label{display:block;margin:14px 0 6px;font-weight:bold}
textarea,input[type=file],input[type=number]{width:100%;padding:10px;border:1px solid var(--line);border-radius:8px;background:#fff;font:inherit}
.buttons{margin-top:16px;display:flex;gap:10px;flex-wrap:wrap}
button,.btn{padding:10px 18px;border:1px solid var(--teal);border-radius:8px;background:#fff;color:var(--navy);font:inherit;font-weight:bold;cursor:pointer;text-decoration:none;display:inline-block}
button.primary{background:var(--teal);color:#fff}summary{cursor:pointer;font-size:18px}
.two{display:grid;grid-template-columns:1fr 1fr;gap:20px}@media(max-width:700px){.two{grid-template-columns:1fr}}
h3.s{color:var(--teal)}h3.f{color:var(--coral)}table{width:100%;border-collapse:collapse;font-size:15px}
th{text-align:left;color:var(--muted);font-weight:normal;vertical-align:top;padding:5px 8px 5px 0;border-top:1px solid var(--line)}
td{padding:5px 8px 5px 0;border-top:1px solid var(--line);vertical-align:top}table.grid th{font-weight:bold;color:var(--navy)}
.chart{width:100%;height:380px}.tag{display:inline-block;background:#e4f4f2;color:#0a5f57;border-radius:6px;padding:1px 8px;margin:2px 3px 2px 0;font-size:14px}
.tag.bad{background:#fdeceb;color:#8a1f17}.tag.b{background:#fbe8e1;color:#8a3318}
@page{margin:12mm}@media print{header,.buttons,.noprint{display:none}details{display:block}details>*{display:block}.chart,tr,.tag,svg{break-inside:avoid;page-break-inside:avoid}h2,h3{break-after:avoid;page-break-after:avoid}.card{overflow:visible}}
</style></head><body>
<header><a href="{{ url_for('index') }}">GeneDisease Finder</a></header><main>
{% with messages = get_flashed_messages() %}{% for m in messages %}<p class="alert">{{ m }}</p>{% endfor %}{% endwith %}
{% if r is none %}
<h1>Enter your genes</h1>
<p class="muted">Paste gene symbols, upload a CSV, TXT or Excel file, or load the example lists.</p>
<form method="post" action="{{ url_for('run') }}" enctype="multipart/form-data">
<label for="genes">Gene set A (commas, spaces or new lines)</label>
<textarea id="genes" name="genes" rows="5" placeholder="TP53, BRCA1, EGFR, KRAS, PTEN"></textarea>
<label for="file">Or upload a file for set A (.csv, .tsv, .txt, .xlsx, max 5 MB)</label>
<input id="file" type="file" name="file" accept=".csv,.tsv,.txt,.xlsx,.xls">
<label for="genes_b">Gene set B (optional, to compare with set A)</label>
<textarea id="genes_b" name="genes_b" rows="3" placeholder="IL6, TNF, INS"></textarea>
<label for="min_score">Minimum disease evidence score (0 to 1)</label>
<input id="min_score" type="number" name="min_score" min="0" max="1" step="0.05" value="0.2" style="max-width:140px">
<label style="font-weight:normal"><input type="checkbox" name="hide_generic" value="1"> Hide very generic diseases (cancer, neoplasm, carcinoma...)</label>
<div class="buttons"><button class="primary" type="submit">Analyze genes</button>
<button type="submit" name="example" value="1">Load example genes</button></div></form>
{% else %}
<h1>Results</h1>
<div class="buttons noprint" style="margin:0 0 12px"><a class="btn" href="{{ url_for('index') }}">New analysis</a>
<button onclick="window.print()">Print or save as PDF</button>
<a class="btn" href="{{ url_for('download', rid=rid, kind='features') }}">Features CSV</a>
<a class="btn" href="{{ url_for('download', rid=rid, kind='pathways') }}">Pathways CSV</a>
<a class="btn" href="{{ url_for('download', rid=rid, kind='diseases') }}">Diseases CSV</a></div>
<div class="card"><p><strong>{{ r.valid_symbols|length }}</strong> valid of <strong>{{ r.input_count }}</strong> genes entered.</p>
{% for n in r.notes %}<p class="muted">{{ n }}</p>{% endfor %}
{% if r.invalid %}<p class="error">Not recognized: {% for i in r.invalid %}{{ i.symbol }}{% if i.suggestion %} (did you mean {{ i.suggestion }}?){% endif %}{% if not loop.last %}, {% endif %}{% endfor %}</p>{% endif %}
{% for w in r.warnings %}<p class="warn">{{ w }}</p>{% endfor %}</div>
{% if r.findings %}<div class="card"><h3 style="margin-top:0">Key findings</h3><ul style="margin:0;padding-left:20px">{% for f in r.findings %}<li>{{ f }}</li>{% endfor %}</ul></div>{% endif %}

<h2>1. Gene features</h2>
{% if r.features_limited %}<p class="muted">Features are shown for the first {{ max_features }} valid genes.</p>{% endif %}
{% for g in r.features %}<details class="card" {% if loop.first %}open{% endif %}>
<summary><strong>{{ g.symbol }}</strong> <span class="muted">{{ g.name }}</span></summary>
{% if g.summary %}<p>{{ g.summary }}</p>{% endif %}
<div class="two"><div><h3 class="s">Structural</h3><table>
<tr><th>Location</th><td>{{ g.structural.map_location or "-" }}</td></tr>
<tr><th>Chromosome</th><td>{{ g.structural.chromosome or "-" }}{% if g.structural.strand %} ({{ g.structural.strand }} strand){% endif %}</td></tr>
<tr><th>Transcripts</th><td>{{ g.structural.transcripts or "-" }}</td></tr>
<tr><th>Max exons</th><td>{{ g.structural.max_exons or "-" }}</td></tr>
<tr><th>Protein length</th><td>{% if g.structural.protein_length %}{{ g.structural.protein_length }} amino acids{% else %}-{% endif %}</td></tr>
<tr><th>Domains</th><td>{% for d in g.structural.domains %}{{ d.name }} ({{ d.start }}-{{ d.end }}){% if not loop.last %}<br>{% endif %}{% else %}-{% endfor %}</td></tr>
<tr><th>ClinVar variants</th><td>{% if g.variants and g.variants.total %}{{ g.variants.total }} in total, {{ g.variants.pathogenic }} pathogenic or likely pathogenic{% else %}-{% endif %}</td></tr>
<tr><th>3D model</th><td>{% if g.structural.alphafold_url %}<a href="{{ g.structural.alphafold_url }}" target="_blank" rel="noopener">AlphaFold DB</a>{% else %}-{% endif %}</td></tr></table></div>
<div><h3 class="f">Functional</h3><table>
<tr><th>Gene type</th><td>{{ g.functional.gene_type or "-" }}</td></tr>
<tr><th>GO process</th><td>{{ g.functional.go_process|join(", ") or "-" }}</td></tr>
<tr><th>GO function</th><td>{{ g.functional.go_function|join(", ") or "-" }}</td></tr>
<tr><th>GO component</th><td>{{ g.functional.go_component|join(", ") or "-" }}</td></tr>
<tr><th>Top tissues (RNA, TPM)</th><td>{% for t in g.expression %}{{ t.tissue }} ({{ t.tpm }}){% if not loop.last %}, {% endif %}{% else %}-{% endfor %}</td></tr>
<tr><th>Druggability</th><td>{% if r.druggable.get(g.symbol) %}{% for l in r.druggable[g.symbol] %}<span class="tag">{{ l }}</span>{% endfor %}{% else %}-{% endif %}</td></tr></table></div></div></details>{% endfor %}

<h2>2. Pathways and common features</h2>
{% if r.pathways %}<div class="card"><div id="c_path" class="chart"></div>
<table class="grid"><tr><th>Pathway</th><th>Genes</th><th>Adj. p-value</th></tr>
{% for p in r.pathways %}<tr><td>{{ p.pathway }}</td><td>{{ p.genes|join(", ") }}</td><td>{{ "%.2e"|format(p.adj_p_value) }}{% if not p.significant %} <span class="muted">(not significant)</span>{% endif %}</td></tr>{% endfor %}</table></div>
{% else %}<p class="muted">No pathways to show.</p>{% endif %}
{% if r.common %}<div class="card"><h3>Common features (shared by 2 or more genes)</h3><table class="grid"><tr><th>Type</th><th>Feature</th><th>Genes</th></tr>
{% for c in r.common %}<tr><td>{{ c.type }}</td><td>{{ c.feature }}</td><td>{{ c.count }}: {{ c.genes|join(", ") }}</td></tr>{% endfor %}</table></div>{% endif %}
{% if r.string and r.string.edges %}<div class="card"><h3>Interaction network (STRING)</h3><div id="c_ppi" class="chart"></div>
<p class="muted">Most connected genes: {% for g, n in r.string.hubs %}{{ g }} ({{ n }}){% if not loop.last %}, {% endif %}{% endfor %}</p></div>{% endif %}

<h2>3. Disease associations</h2>
{% if r.ranked %}<div class="card"><div id="c_dis" class="chart"></div>
<form method="get" class="noprint" action="{{ url_for('show', rid=rid) }}">
<label for="ms">Evidence score filter: minimum score <strong id="msv">{{ "%.2f"|format(r.min_score) }}</strong> (applies to the tables and charts below)</label>
<input id="ms" type="range" name="min_score" min="0" max="0.9" step="0.05" value="{{ r.min_score }}" oninput="document.getElementById('msv').textContent=parseFloat(this.value).toFixed(2)" onchange="this.form.submit()" style="width:100%;max-width:420px">
<label style="font-weight:normal"><input type="checkbox" name="hide_generic" value="1" {% if r.hide_generic %}checked{% endif %} onchange="this.form.submit()"> Hide very generic diseases (cancer, neoplasm, carcinoma...)</label></form>
<table class="grid"><tr><th>#</th><th>Disease</th><th>Genes</th><th>Mean score</th></tr>
{% for d in r.ranked[:15] %}<tr><td>{{ loop.index }}</td><td>{{ d.disease }}</td><td>{{ d.n_genes }}: {{ d.genes|join(", ") }}</td><td>{{ "%.2f"|format(d.mean_score) }}</td></tr>{% endfor %}</table></div>
<div class="card"><h3>Gene-disease heatmap (association score)</h3><div id="c_heat" class="chart" style="height:460px"></div></div>
<div class="card"><h3>Gene-disease network</h3><div id="c_gd" class="chart" style="height:460px"></div></div>
{% if r.related %}<div class="card"><h3>Related diseases (genes shared with {{ r.ranked[0].disease }})</h3>
{% for x in r.related %}<span class="tag b">{{ x.disease }} ({{ x.shared }})</span>{% endfor %}</div>{% endif %}
{% if r.druggable %}<div class="card"><h3>Drug targets</h3><table class="grid"><tr><th>Gene</th><th>Druggability evidence</th></tr>
{% for g, labels in r.druggable.items() %}<tr><td>{{ g }}</td><td>{{ labels|join(", ") }}</td></tr>{% endfor %}</table></div>{% endif %}
{% else %}<p class="muted">No disease associations to show.</p>{% endif %}

{% if r.compare %}<h2>4. Compare gene sets</h2><div class="card">
<p>Set A: <strong>{{ r.compare.count_a }}</strong> genes. Set B: <strong>{{ r.compare.count_b }}</strong> genes.
{% if r.compare.invalid_b %}<span class="error">Not recognized in B: {{ r.compare.invalid_b|join(", ") }}</span>{% endif %}</p>
<svg viewBox="0 0 600 240" width="100%" style="max-width:600px;display:block;margin:0 auto" role="img" aria-label="Genes only in set A, shared, and only in set B">
<circle cx="240" cy="110" r="95" fill="rgba(15,148,136,0.3)" stroke="#0f9488"/><circle cx="360" cy="110" r="95" fill="rgba(212,85,47,0.3)" stroke="#d4552f"/>
<g text-anchor="middle" fill="#0b2545" font-family="Calibri, Arial">
<text x="200" y="108" font-size="28" font-weight="bold">{{ r.compare.only_a|length }}</text><text x="200" y="130" font-size="14">only A</text>
<text x="300" y="108" font-size="28" font-weight="bold">{{ r.compare.shared|length }}</text><text x="300" y="130" font-size="14">shared</text>
<text x="400" y="108" font-size="28" font-weight="bold">{{ r.compare.only_b|length }}</text><text x="400" y="130" font-size="14">only B</text>
<text x="200" y="228" font-size="15" font-weight="bold">Set A</text><text x="400" y="228" font-size="15" font-weight="bold">Set B</text></g></svg>
<table class="grid"><tr><th>Shared genes</th><td>{{ r.compare.shared|join(", ") or "none" }}</td></tr>
<tr><th>Only in A</th><td>{{ r.compare.only_a|join(", ") or "none" }}</td></tr>
<tr><th>Only in B</th><td>{{ r.compare.only_b|join(", ") or "none" }}</td></tr>
{% if r.compare.shared_pathways is defined %}<tr><th>Shared pathways</th><td>{{ r.compare.shared_pathways|join(", ") or "none" }}</td></tr>{% endif %}</table></div>{% endif %}

<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script>
const D = {{ charts|tojson }};
const cfg = {displayModeBar:false, responsive:true};
const lay = (extra) => Object.assign({margin:{l:10,r:20,t:30,b:40}, paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)', font:{family:'Calibri, Arial'}}, extra||{});
function bar(id, items, color, title, xtitle){
  const el=document.getElementById(id); if(!el||!items||!items.length) return;
  Plotly.newPlot(el,[{type:'bar',orientation:'h',x:items.map(i=>i.value).reverse(),y:items.map(i=>i.label).reverse(),marker:{color:color},text:items.map(i=>i.value).reverse(),textposition:'outside'}],
    lay({title:{text:title,font:{size:15}},margin:{l:260,r:30,t:40,b:30},xaxis:{title:xtitle,dtick:1}}),cfg);
}
bar('c_path', D.pathways, '#0f9488', 'Input genes in each enriched pathway', 'Number of input genes in the pathway');
bar('c_dis', D.diseases, '#d4552f', 'Genes supporting each disease', 'Number of input genes linked to the disease');
if (D.heatmap && document.getElementById('c_heat')) {
  Plotly.newPlot('c_heat',[{type:'heatmap',x:D.heatmap.genes,y:D.heatmap.diseases,z:D.heatmap.z,colorscale:[[0,'#e8f4f2'],[1,'#0b6b62']],zmin:0,zmax:1,hoverongaps:false,colorbar:{title:'Score'}}],
    lay({margin:{l:260,r:20,t:20,b:80},xaxis:{tickangle:-45}}),cfg);
}
if (D.gene_disease && document.getElementById('c_gd')) {
  const g=D.gene_disease, gy={}, dy={};
  g.genes.forEach((n,i)=>gy[n]=-i/Math.max(g.genes.length-1,1)); g.diseases.forEach((n,i)=>dy[n]=-i/Math.max(g.diseases.length-1,1));
  const lx=[], ly=[]; g.edges.forEach(e=>{lx.push(0,1,null); ly.push(gy[e[0]],dy[e[1]],null);});
  Plotly.newPlot('c_gd',[
    {type:'scatter',mode:'lines',x:lx,y:ly,line:{color:'#b8c4cb',width:1},hoverinfo:'skip'},
    {type:'scatter',mode:'markers+text',x:g.genes.map(()=>0),y:g.genes.map(n=>gy[n]),text:g.genes,textposition:'middle left',cliponaxis:false,marker:{size:12,color:'#0f9488'},hoverinfo:'text'},
    {type:'scatter',mode:'markers+text',x:g.diseases.map(()=>1),y:g.diseases.map(n=>dy[n]),text:g.diseases,textposition:'middle right',cliponaxis:false,marker:{size:14,color:'#d4552f'},hoverinfo:'text'}],
    lay({showlegend:false,margin:{l:90,r:260,t:10,b:10},xaxis:{visible:false,range:[-0.2,1.2]},yaxis:{visible:false,range:[-1.1,0.1]}}),cfg);
}
if (D.ppi && document.getElementById('c_ppi')) {
  const p=D.ppi, ex=[], ey=[]; p.edges.forEach(e=>{ex.push(p.pos[e[0]][0],p.pos[e[1]][0],null); ey.push(p.pos[e[0]][1],p.pos[e[1]][1],null);});
  const nodes=Object.keys(p.pos);
  Plotly.newPlot('c_ppi',[
    {type:'scatter',mode:'lines',x:ex,y:ey,line:{color:'#b8c4cb',width:1},hoverinfo:'skip'},
    {type:'scatter',mode:'markers+text',x:nodes.map(n=>p.pos[n][0]),y:nodes.map(n=>p.pos[n][1]),text:nodes,textposition:'top center',cliponaxis:false,marker:{size:12,color:'#0f9488'},hoverinfo:'text'}],
    lay({showlegend:false,margin:{l:20,r:20,t:20,b:20},xaxis:{visible:false},yaxis:{visible:false,scaleanchor:'x'}}),cfg);
}
</script>
{% endif %}
</main></body></html>"""


def render_page(result=None, rid=None):
    return render_template_string(
        PAGE, r=result, rid=rid, max_features=MAX_FEATURE_GENES,
        charts=(result or {}).get("charts") or {})


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
