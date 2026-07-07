import pandas as pd

df = pd.read_csv(r'C:\Users\seck-s\Documents\_UNC\Projets\Article\RESULTAT\Figure5_spectral_candidates.csv')
print('Shape:', df.shape)
print()
print('category_spectral:')
print(df['category_spectral'].value_counts())
print()
print('found_mass:')
print(df['found_mass'].value_counts())
print()
ms = df[df['found_mass'].str.lower()=='yes']
dup = ms['exact_mass'].value_counts()
print('Groupes isomeres (meme masse, MS confirmes):', (dup>1).sum(), 'masses,', dup[dup>1].sum(), 'molecules')
print(dup[dup>1].head(10))
print()
print('combined_score stats:')
print(df['combined_score'].describe())
