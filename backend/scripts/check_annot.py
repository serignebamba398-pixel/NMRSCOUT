import pandas as pd, csv

def read_auto(path):
    with open(path, 'r', encoding='utf-8-sig', errors='ignore') as f:
        sample = f.read(2048)
    try:
        sep = csv.Sniffer().sniff(sample, delimiters=[',',';','\t']).delimiter
    except:
        sep = ','
    return pd.read_csv(path, sep=sep, engine='python')

BASE = r"C:\Users\seck-s\Documents\_UNC\Projets\Article\RESULTAT"

sp = read_auto(f"{BASE}/Figure5_spectral_candidates.csv")
h1 = read_auto(f"{BASE}/annotation_1H_latest.csv")

print("spectral Source sample:", sp['Source'].iloc[:5].tolist())
print("1H Source sample:", h1['Source'].iloc[:5].tolist())
print()

# Voir si les noms matchent directement ou avec suffix _predicted
sp_names = set(sp['Source'])
h1_names = set(h1['Source'])

direct = sp_names & h1_names
predicted = {n+'_predicted' for n in sp_names} & h1_names

print(f"Match direct Source→Source : {len(direct)}")
print(f"Match avec suffix _predicted : {len(predicted)}")
print()
print("Colonnes 1H :", list(h1.columns[:8]))
