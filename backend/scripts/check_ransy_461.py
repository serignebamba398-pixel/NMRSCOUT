import os, csv
import pandas as pd

FOLDER = r"C:\Users\seck-s\Documents\_UNC\Projets\Article\RESULTAT"

def read_auto(path):
    with open(path, 'r', encoding='utf-8-sig', errors='ignore') as f:
        sample = f.read(2048)
    try:
        sep = csv.Sniffer().sniff(sample, delimiters=[',', ';', '\t']).delimiter
    except:
        sep = ','
    return pd.read_csv(path, sep=sep, engine='python')

mols = ['DNP000042', 'NPA007877']

for method, fname in [('1H','annotation_1H_latest.csv'),('13C','annotation_13C_latest.csv'),
                       ('HSQC','annotation_HSQC_latest.csv'),('COSY','annotation_COSY_latest.csv')]:
    df = read_auto(os.path.join(FOLDER, fname))
    print(f"\n{'='*60}")
    print(f"MÉTHODE : {method}")
    print(f"{'='*60}")
    for mol in mols:
        row = df[df['Source'] == mol + '_predicted']
        if len(row) == 0:
            print(f"\n  {mol} : NON TROUVÉ")
            continue
        row = row.iloc[0]
        print(f"\n  {mol} (score={row.get('intersection_coverage' if method in ('1H','13C') else 'coverage_ransy', 'N/A')}) :")
        print(f"    ppm_molecule          : {row.get('ppm_molecule','')}")
        print(f"    RANSY_signals         : {row.get('RANSY_signals','')}")
        print(f"    ppm_molecule_matched  : {row.get('ppm_molecule_matched','')}")
        print(f"    intersection_ppms     : {row.get('intersection_ppms','')}")
        print(f"    anchor_shifts         : {row.get('anchor_shifts','')}")
