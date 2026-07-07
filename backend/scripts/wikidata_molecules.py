#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wikidata_molecules.py  (version corrigée)
-----------------------------------------
Récupère depuis Wikidata des molécules et leurs informations
(nom, SMILES, InChIKey, formule, PubChem CID) ainsi que l'organisme
source (espèce, genre, famille) via la propriété P703 (« présent dans le taxon »).

S'appuie sur les données LOTUS intégrées à Wikidata (produits naturels).

Quatre modes :
    --inchikey BLUAFEHZUWYNDE-NNWCWBAJSA-N      une molécule précise
    --species  "Centella asiatica"              toutes les molécules d'une espèce
    --genus    "Centella"                       toutes les molécules d'un genre (résout le QID automatiquement)
    --clade    Q25400                           toutes les molécules sous un clade

Sortie : CSV (.csv) ou Excel (.xlsx) selon l'extension de -o.

Dépendances :
    pip install requests pandas
    pip install openpyxl        # uniquement pour la sortie .xlsx

Exemples :
    python wikidata_molecules.py --species "Centella asiatica" -o centella.csv
    python wikidata_molecules.py --genus "Centella" -o centella_genre.csv
    python wikidata_molecules.py --species "Centella asiatica" -o centella.xlsx
    python wikidata_molecules.py --clade Q25400 --page-size 1000 -o asteraceae.csv
"""

import os
import sys
import time
import argparse
from typing import List, Dict, Optional

import requests
import pandas as pd

WIKIDATA_SPARQL = "https://query.wikidata.org/sparql"

# User-Agent explicite : Wikidata l'exige. Remplacez l'email par le vôtre.
HEADERS = {
    "User-Agent": "MoleculeTaxonomyFetcher/1.1 (research; contact@example.org)",
    "Accept": "application/sparql-results+json",
}


# ===================================================================
# 1a) RÉSOLUTION AUTOMATIQUE DU QID D'UN GENRE
# ===================================================================
def resolve_genus_qid(genus_name: str) -> str:
    """Trouve le QID Wikidata d'un genre taxonomique par son nom.

    Cherche un taxon de rang 'genre' (Q34740) ayant le nom scientifique donné,
    et qui a au moins une molécule liée (P703) pour éviter les doublons vides.
    """
    safe = genus_name.replace('"', '\\"')
    query = f"""SELECT ?genus (COUNT(DISTINCT ?compound) AS ?n) WHERE {{
  ?genus wdt:P105 wd:Q34740 .
  ?genus wdt:P225 "{safe}" .
  ?taxon wdt:P171* ?genus .
  ?compound wdt:P703 ?taxon .
}}
GROUP BY ?genus
ORDER BY DESC(?n)
LIMIT 5
"""
    print(f"Recherche du genre '{genus_name}' sur Wikidata...", file=sys.stderr)
    bindings = run_query(query)
    if not bindings:
        print(f"  [ERREUR] Aucun genre '{genus_name}' trouve sur Wikidata.", file=sys.stderr)
        sys.exit(1)

    # Prendre celui avec le plus de molécules
    best = bindings[0]
    qid = best["genus"]["value"].rsplit("/", 1)[-1]
    n_molecules = int(best["n"]["value"])
    print(f"  Genre '{genus_name}' -> {qid} ({n_molecules} molecules)", file=sys.stderr)

    if len(bindings) > 1:
        print(f"  [INFO] {len(bindings)} QIDs trouves pour '{genus_name}', "
              f"utilisation de {qid} (le plus de molecules).", file=sys.stderr)
        for b in bindings:
            q = b["genus"]["value"].rsplit("/", 1)[-1]
            n = b["n"]["value"]
            print(f"    - {q} : {n} molecules", file=sys.stderr)

    return qid


# ===================================================================
# 1b) CONSTRUCTION DES REQUÊTES SPARQL  (corrigée : plus de BIND invalide)
# ===================================================================
def build_query(mode: str, value: str, limit: int, offset: int) -> str:
    # corps commun : on utilise partout la même variable ?taxon
    common = """
  ?taxon wdt:P225 ?taxonName .
  OPTIONAL { ?compound wdt:P233 ?smiles. }
  OPTIONAL { ?compound wdt:P235 ?inchikey. }
  OPTIONAL { ?compound wdt:P274 ?formula. }
  OPTIONAL { ?compound wdt:P662 ?pubchem. }
  OPTIONAL { ?taxon wdt:P171* ?g. ?g wdt:P105 wd:Q34740. ?g wdt:P225 ?genus. }
  OPTIONAL { ?taxon wdt:P171* ?f. ?f wdt:P105 wd:Q35409. ?f wdt:P225 ?family. }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "fr,en". }
"""

    safe = value.replace('"', '\\"')  # échappe d'éventuels guillemets

    if mode == "inchikey":
        head = (f'  ?compound wdt:P235 "{safe}" .\n'
                f'  ?compound wdt:P703 ?taxon .')
    elif mode == "species":
        head = (f'  ?taxon wdt:P225 "{safe}" .\n'
                f'  ?compound wdt:P703 ?taxon .')
    elif mode == "clade":
        head = (f'  VALUES ?clade {{ wd:{value} }}\n'
                f'  ?taxon wdt:P171* ?clade .\n'
                f'  ?compound wdt:P703 ?taxon .')
    else:
        raise ValueError(f"Mode inconnu : {mode}")

    return f"""SELECT DISTINCT ?compound ?compoundLabel ?smiles ?inchikey ?formula ?pubchem ?taxonName ?genus ?family
WHERE {{
{head}
{common}
}}
LIMIT {limit}
OFFSET {offset}
"""


# ===================================================================
# 2) EXÉCUTION D'UNE REQUÊTE (retry sur timeout / surcharge)
# ===================================================================
def run_query(query: str, max_retries: int = 3, backoff: float = 5.0) -> Optional[List[Dict]]:
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(
                WIKIDATA_SPARQL,
                params={"query": query, "format": "json"},
                headers=HEADERS,
                timeout=90,
            )
            if resp.status_code == 200:
                return resp.json()["results"]["bindings"]
            if resp.status_code in (429, 500, 502, 503, 504):
                wait = backoff * attempt
                print(f"  [HTTP {resp.status_code}] tentative {attempt}/{max_retries}, "
                      f"attente {wait:.0f}s...", file=sys.stderr)
                if attempt == 1:
                    print(f"  detail: {resp.text[:200]}", file=sys.stderr)
                time.sleep(wait)
                continue
            print(f"  [ERREUR HTTP {resp.status_code}] {resp.text[:300]}", file=sys.stderr)
            return None
        except requests.exceptions.Timeout:
            wait = backoff * attempt
            print(f"  [TIMEOUT] tentative {attempt}/{max_retries}, attente {wait:.0f}s...",
                  file=sys.stderr)
            time.sleep(wait)
        except requests.exceptions.RequestException as e:
            print(f"  [EXCEPTION] {e}", file=sys.stderr)
            time.sleep(backoff)
    print("  [ÉCHEC] abandon après plusieurs tentatives.", file=sys.stderr)
    return None


# ===================================================================
# 3) APLATISSAGE
# ===================================================================
def flatten(bindings: List[Dict]) -> List[Dict]:
    rows = []
    for b in bindings:
        rows.append({
            "qid":        b.get("compound", {}).get("value", "").rsplit("/", 1)[-1],
            "nom":        b.get("compoundLabel", {}).get("value", ""),
            "SMILES":     b.get("smiles", {}).get("value", ""),
            "InChIKey":   b.get("inchikey", {}).get("value", ""),
            "formule":    b.get("formula", {}).get("value", ""),
            "PubChemCID": b.get("pubchem", {}).get("value", ""),
            "espece":     b.get("taxonName", {}).get("value", ""),
            "genre":      b.get("genus", {}).get("value", ""),
            "famille":    b.get("family", {}).get("value", ""),
        })
    return rows


# ===================================================================
# 4) RÉCUPÉRATION COMPLÈTE (pagination)
# ===================================================================
def fetch_all(mode: str, value: str, page_size: int = 1000,
              max_pages: int = 50, pause: float = 1.0) -> pd.DataFrame:
    all_rows: List[Dict] = []
    offset = 0
    pages = 1 if mode == "inchikey" else max_pages

    for page in range(pages):
        print(f"[Page {page + 1}] offset={offset} ...", file=sys.stderr)
        bindings = run_query(build_query(mode, value, page_size, offset))
        if bindings is None:
            print("  Arrêt : échec de la requête.", file=sys.stderr)
            break
        if not bindings:
            print("  Plus de résultats.", file=sys.stderr)
            break

        rows = flatten(bindings)
        all_rows.extend(rows)
        print(f"  +{len(rows)} lignes (total {len(all_rows)})", file=sys.stderr)

        if len(bindings) < page_size:
            break
        offset += page_size
        time.sleep(pause)

    df = pd.DataFrame(all_rows)
    if not df.empty:
        df.drop_duplicates(inplace=True)
        df.sort_values(by=["nom", "espece"], inplace=True, ignore_index=True)
    return df


# ===================================================================
# 5) ÉCRITURE (CSV ou XLSX selon l'extension)
# ===================================================================
def write_output(df: pd.DataFrame, output: str):
    ext = os.path.splitext(output)[1].lower()
    if ext in (".xlsx", ".xls"):
        try:
            df.to_excel(output, index=False, engine="openpyxl")
        except ModuleNotFoundError:
            print("  [INFO] openpyxl manquant -> installez-le avec : pip install openpyxl",
                  file=sys.stderr)
            fallback = os.path.splitext(output)[0] + ".csv"
            df.to_csv(fallback, index=False, encoding="utf-8-sig")
            print(f"  Sortie repliée en CSV : {fallback}", file=sys.stderr)
            return fallback
    else:
        df.to_csv(output, index=False, encoding="utf-8-sig")
    return output


# ===================================================================
# 6) CLI
# ===================================================================
def main():
    p = argparse.ArgumentParser(
        description="Récupère molécules + taxonomie source depuis Wikidata (LOTUS)."
    )
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--inchikey", help="InChIKey d'une molécule précise")
    g.add_argument("--species",  help="Nom scientifique d'une espèce (ex. 'Centella asiatica')")
    g.add_argument("--genus",    help="Nom d'un genre (ex. 'Centella') - resout le QID automatiquement")
    g.add_argument("--clade",    help="QID d'un clade (ex. Q25400 pour Asteraceae)")

    p.add_argument("-o", "--output", default="molecules_wikidata.csv",
                   help="Fichier de sortie (.csv ou .xlsx)")
    p.add_argument("--page-size", type=int, default=1000,
                   help="Taille de page (LIMIT). Réduisez en cas de timeout.")
    p.add_argument("--max-pages", type=int, default=50,
                   help="Nombre maximum de pages.")
    args = p.parse_args()

    if args.inchikey:
        mode, value = "inchikey", args.inchikey
    elif args.species:
        mode, value = "species", args.species
    elif args.genus:
        qid = resolve_genus_qid(args.genus)
        mode, value = "clade", qid
    else:
        mode, value = "clade", args.clade

    print(f"Mode : {mode} = {value}", file=sys.stderr)
    df = fetch_all(mode, value, page_size=args.page_size, max_pages=args.max_pages)

    if df.empty:
        print("Aucun résultat.", file=sys.stderr)
        sys.exit(1)

    out = write_output(df, args.output)
    print(f"\n✓ {len(df)} lignes écrites dans {out}", file=sys.stderr)
    with pd.option_context("display.max_columns", None, "display.width", 160):
        print(df.head(10).to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()