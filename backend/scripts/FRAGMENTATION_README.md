# Fragmentation Hamiltonienne pour Systèmes de Spins

## Vue d'ensemble

Implémentation de l'algorithme de fragmentation de Castillo et al. (2011) pour simuler de grands systèmes de spins (> 10 spins) en RMN.

**Référence**: Castillo, Patiny & Wist, *J. Magn. Reson.* **209**, 123-130 (2011)
"Fast and Accurate Algorithm for the Simulation of NMR spectra of Large Spin Systems"

## Principe

Pour un système de N spins couplés, la simulation quantique exacte nécessite la diagonalisation d'une matrice Hamiltonienne de taille 2^N × 2^N:

| N spins | Dimension | Éléments | Faisabilité |
|---------|-----------|----------|-------------|
| 5       | 32×32     | 1,024    | ✓ Rapide    |
| 10      | 1024×1024 | 1M       | ✓ Limite    |
| 15      | 32768×... | 1G       | ✗ Impossible|
| 20      | 1M×...    | 1T       | ✗ Impossible|

**Solution**: Fragmenter le système en sous-systèmes faiblement couplés.

### Critère de fragmentation

Deux spins i et j sont **faiblement couplés** (peuvent être séparés) si:

```
|J_ij / Δν_ij| < seuil
```

Où:
- `J_ij` = constante de couplage scalaire (Hz)
- `Δν_ij` = différence de fréquence de Larmor (Hz)
- `seuil` = typiquement 0.05

Si cette condition est vérifiée, les deux spins peuvent être simulés dans des fragments séparés.

## Algorithme

### Étape 1: Construction du graphe de couplage

```python
def _build_coupling_graph(frequencies_hz, J_matrix, threshold_ratio=0.05):
    """
    Construit un graphe où:
    - Nœuds = spins
    - Arête (i,j) existe si |J_ij/Δν_ij| > threshold_ratio
    """
```

### Étape 2: Détection des fragments

```python
def _find_spin_fragments(frequencies_hz, J_matrix, threshold_ratio=0.05):
    """
    Retourne les composantes connexes du graphe.
    Ex: [[0,1,2], [3,4,5,6], [7]] = 3 fragments
    """
```

### Étape 3: Simulation de chaque fragment

```python
def _simulate_fragment(fragment_indices, frequencies_hz, J_matrix):
    """
    Simule un fragment isolé par diagonalisation exacte.
    Limite: ≤ 10 spins par fragment
    """
```

### Étape 4: Combinaison des spectres

```python
def _simulate_large_spin_system(frequencies_hz, J_matrix, threshold_ratio=0.05):
    """
    Pipeline complet:
    1. Détecter les fragments
    2. Simuler chaque fragment
    3. Combiner les transitions
    """
```

## Intégration

La fragmentation est intégrée directement dans `_simulate_spin_system_exact()`:

```python
def _simulate_spin_system_exact(
    frequencies_hz,
    J_matrix,
    intensity_threshold=0.001,
    use_fragmentation=True,      # ← Nouveau paramètre
    fragmentation_threshold=0.05  # ← Seuil |J/Δν|
):
    """
    Simule un système de spins:
    - ≤ 10 spins: simulation directe
    - > 10 spins + use_fragmentation=True: tentative de fragmentation
    - Échec: fallback vers approximation premier ordre
    """
```

## Configuration

Dans la classe `SimConfig`:

```python
@dataclass
class SimConfig:
    # ... autres options ...

    # Fragmentation de systèmes de spins (Castillo et al. 2011)
    use_spin_fragmentation: bool = True          # Activer la fragmentation
    fragmentation_threshold_ratio: float = 0.05  # Seuil |J/Δν|
```

## Exemples d'utilisation

### Cas 1: Système fragmentable (2 paires AB éloignées)

```python
# 4 spins: 2 systèmes AB indépendants
frequencies_hz = [
    800,  1000,  # Paire 1: δ ≈ 2.0-2.5 ppm @ 400 MHz
    2800, 3000   # Paire 2: δ ≈ 7.0-7.5 ppm @ 400 MHz
]

J_matrix = np.array([
    [0, 7, 0, 0],  # Couplage intra-paire
    [7, 0, 0, 0],
    [0, 0, 0, 7],
    [0, 0, 7, 0]
])

# Fragmentation automatique: 2 fragments de 2 spins
transitions = _simulate_spin_system_exact(
    frequencies_hz, J_matrix,
    use_fragmentation=True
)
```

**Résultat**: 2 fragments détectés → simulation réussie

### Cas 2: Système non fragmentable (tous couplés)

```python
# 12 spins tous proches et couplés
frequencies_hz = [1000 + i*5 for i in range(12)]  # Très proches

J_matrix = np.ones((12, 12)) * 7.0  # Tous couplés

# Fragmentation échouera → fallback approximation 1er ordre
transitions = _simulate_spin_system_exact(
    frequencies_hz, J_matrix,
    use_fragmentation=True
)
```

**Résultat**: 1 fragment de 12 spins (trop grand) → fallback vers approximation

## Cas d'usage réels

### ✓ Fonctionne bien pour:

1. **Molécules avec groupes chimiquement distincts**
   - Exemple: Composé aromatique + chaîne alkyle
   - Fragment 1: Aromatique (δ ≈ 7 ppm)
   - Fragment 2: Alkyle (δ ≈ 1-2 ppm)

2. **Systèmes hétérocycliques**
   - Cycles indépendants avec peu d'interactions

3. **Stéroïdes avec substituants**
   - Corps stéroïde + chaîne latérale

### ✗ Limitations:

1. **Longues chaînes alkyle homogènes**
   - Exemple: Octane CH₃-(CH₂)₆-CH₃
   - Problème: CH₂ centraux quasi-dégénérés (δ ≈ 1.26 ppm)
   - Résultat: Un seul fragment > 10 spins → non fragmentable

2. **Systèmes fortement couplés étendus**
   - Cycles fusionnés avec couplages longue distance
   - Systèmes aromatiques polycycliques

3. **Protons quasi-dégénérés en grand nombre**
   - Δν ≈ 0 → Impossible de fragmenter

**Solution pour ces cas**: Spinach (RSS) ou approximation premier ordre

## Dépendances

### networkx (optionnel)

La fragmentation utilise `networkx` pour l'analyse de graphe (composantes connexes).

**Si networkx n'est pas installé:**
- La fragmentation est automatiquement désactivée
- Fallback vers:
  - Simulation directe si ≤ 10 spins
  - Approximation premier ordre si > 10 spins

**Installation** (si nécessaire):
```bash
pip install networkx
```

## Complexité et performance

### Sans fragmentation:
- N spins → matrice 2^N × 2^N
- Temps: O(8^N) pour la diagonalisation
- Limite pratique: N ≤ 10

### Avec fragmentation:
- F fragments de taille N_i (i=1..F)
- Temps: Σ O(8^N_i)
- Gain: Si 12 spins → 2 fragments de 6 → 64×64 × 2 au lieu de 4096×4096

**Exemple**:
```
Sans fragmentation: 12 spins → 4096×4096 = 16M éléments
Avec fragmentation: 2×6 spins → 2×(64×64) = 8K éléments
Gain: ~2000×
```

## Tests

Exécuter la démonstration conceptuelle:
```bash
cd backend/scripts
python test_fragmentation_example.py
```

## Logs et débogage

Le code affiche des informations dans stderr:

```
[INFO] Système de 12 spins > 10, tentative de fragmentation
[INFO] Fragmentation: 12 spins, seuil J/Δν = 0.05
[INFO] → 3 fragments détectés
[INFO]   Fragment 1: 4 spins [0, 1, 2, 3]
[INFO]   Fragment 2: 5 spins [4, 5, 6, 7, 8]
[INFO]   Fragment 3: 3 spins [9, 10, 11]
[INFO]   Fragment 1: 24 transitions
[INFO]   Fragment 2: 48 transitions
[INFO]   Fragment 3: 12 transitions
[INFO] Total: 84 transitions combinées
```

En cas d'échec:
```
[WARN] Fragment maximum de 12 spins > 10, fragmentation insuffisante
[WARN] → Augmenter threshold_ratio ou utiliser approximation premier ordre
[WARN] Fragmentation échouée, fallback approximation premier ordre
```

## Références

1. **Castillo, A. M., Patiny, L., & Wist, J.** (2011). Fast and accurate algorithm for the simulation of NMR spectra of large spin systems. *Journal of Magnetic Resonance*, 209(2), 123-130.

2. **Pople, J. A., Schneider, W. G., & Bernstein, H. J.** (1959). *High-resolution nuclear magnetic resonance*. McGraw-Hill.

3. **Documentation hamiltonien.md** - Détails de l'implémentation Python

## Auteur

Implémentation: Décembre 2024
Basée sur la documentation hamiltonien.md
