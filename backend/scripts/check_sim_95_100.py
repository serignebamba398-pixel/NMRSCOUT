import os, csv
import numpy as np
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

def parse_1d(s):
    if pd.isna(s) or not str(s).strip(): return []
    return [float(x) for x in str(s).split('|') if x.strip()]

def parse_2d(s):
    if pd.isna(s) or not str(s).strip(): return []
    out = []
    for pair in str(s).split('||'):
        parts = pair.split('|')
        if len(parts) == 2:
            try: out.append((float(parts[0]), float(parts[1])))
            except: pass
    return list(set(out))

def tanimoto_1d(A, B, tol):
    if not A or not B: return 0.0
    A_, B_ = sorted(set(A)), sorted(set(B))
    matched, Bc = 0, B_.copy()
    for a in A_:
        diffs = [abs(b - a) for b in Bc]
        if not diffs: break
        i = min(range(len(diffs)), key=lambda i: diffs[i])
        if diffs[i] <= tol:
            matched += 1; Bc.pop(i)
    union = len(set(A_)) + len(set(B_)) - matched
    return matched / union if union > 0 else 0.0

def tanimoto_2d(A, B, t1, t2):
    if not A or not B: return 0.0
    matched, used = 0, set()
    for a in A:
        for i, b in enumerate(B):
            if i not in used and abs(b[0]-a[0])<=t1 and abs(b[1]-a[1])<=t2:
                matched += 1; used.add(i); break
    union = len(A) + len(B) - matched
    return matched / union if union > 0 else 0.0

TOL = {'1H': 0.15, '13C': 2.0, 'HSQC': (3.0, 0.2), 'COSY': (0.2, 0.2)}

# Charger shifts
shifts = {}
for method, fname in [('1H','annotation_1H_latest.csv'),('13C','annotation_13C_latest.csv'),
                       ('HSQC','annotation_HSQC_latest.csv'),('COSY','annotation_COSY_latest.csv')]:
    df = read_auto(os.path.join(FOLDER, fname))
    for _, row in df.iterrows():
        mol = str(row['Source']).replace('_predicted','')
        if mol not in shifts:
            shifts[mol] = {'1H':[],'13C':[],'HSQC':[],'COSY':[]}
        if method in ('1H','13C'):
            shifts[mol][method] = parse_1d(row.get('ppm_molecule',''))
        else:
            shifts[mol][method] = parse_2d(row.get('ppm_molecule',''))

# Charger le CSV smart isomers
df = read_auto(os.path.join(FOLDER, 'Figure5_smart_isomers_candidates.csv'))

# Filtrer best + close dans le bin 95-100
high = df[(df['combined_score'] >= 95) &
          (df['nmr_iso_category'].isin(['MS_best_isomer', 'MS_close_isomer', 'MS_unique']))]

print(f"Molécules MS avec score 95-100 : {len(high)}")
print(high[['Source','combined_score','nmr_iso_category','exact_mass']].to_string())
print()

# Pour chaque groupe de masse dans ce bin, calculer la sim entre membres
for mass_val, group in high.groupby('exact_mass'):
    if len(group) < 2:
        continue
    mols = group[['Source','combined_score','nmr_iso_category']].values
    print(f"\nMasse {mass_val} — {len(group)} molécules :")
    for i in range(len(mols)):
        for j in range(i+1, len(mols)):
            m1, s1, c1 = mols[i]
            m2, s2, c2 = mols[j]
            sh1 = shifts.get(m1, {'1H':[],'13C':[],'HSQC':[],'COSY':[]})
            sh2 = shifts.get(m2, {'1H':[],'13C':[],'HSQC':[],'COSY':[]})
            t1 = tanimoto_1d(sh1['1H'],  sh2['1H'],  TOL['1H'])
            t2 = tanimoto_1d(sh1['13C'], sh2['13C'], TOL['13C'])
            t3 = tanimoto_2d(sh1['HSQC'],sh2['HSQC'],TOL['HSQC'][0],TOL['HSQC'][1])
            t4 = tanimoto_2d(sh1['COSY'],sh2['COSY'],TOL['COSY'][0],TOL['COSY'][1])
            mean = (t1+t2+t3+t4)/4
            print(f"  {m1}({c1}, {s1:.1f}) vs {m2}({c2}, {s2:.1f})")
            print(f"    1H={t1*100:.1f}%  13C={t2*100:.1f}%  HSQC={t3*100:.1f}%  COSY={t4*100:.1f}%  → moy={mean*100:.1f}%")

# Stats globales pour tous les groupes avec best+close
print("\n\n=== STATS GLOBALES : similarité best vs close isomers ===")
sims = []
for mass_val, group in df[df['nmr_iso_category'].isin(['MS_best_isomer','MS_close_isomer'])].groupby('exact_mass'):
    if len(group) < 2: continue
    best_row = group[group['nmr_iso_category']=='MS_best_isomer']
    close_rows = group[group['nmr_iso_category']=='MS_close_isomer']
    if len(best_row)==0: continue
    best_mol = best_row.iloc[0]['Source']
    sh_best = shifts.get(best_mol, {'1H':[],'13C':[],'HSQC':[],'COSY':[]})
    for _, cr in close_rows.iterrows():
        sh_c = shifts.get(cr['Source'], {'1H':[],'13C':[],'HSQC':[],'COSY':[]})
        t1 = tanimoto_1d(sh_c['1H'],  sh_best['1H'],  TOL['1H'])
        t2 = tanimoto_1d(sh_c['13C'], sh_best['13C'], TOL['13C'])
        t3 = tanimoto_2d(sh_c['HSQC'],sh_best['HSQC'],TOL['HSQC'][0],TOL['HSQC'][1])
        t4 = tanimoto_2d(sh_c['COSY'],sh_best['COSY'],TOL['COSY'][0],TOL['COSY'][1])
        sims.append((t1+t2+t3+t4)/4*100)

if sims:
    sims = np.array(sims)
    print(f"N paires best-close : {len(sims)}")
    print(f"Similarité spectrale moyenne : {sims.mean():.1f}%")
    print(f"Médiane : {np.median(sims):.1f}%")
    print(f"Min/Max : {sims.min():.1f}% / {sims.max():.1f}%")
