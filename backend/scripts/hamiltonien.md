La fragmentation suit l'algorithme de Castillo et al. (2011) en 4 étapes. Voici comment l'implémenter directement dans votre code Python :

## Principe

```
Graphe de couplage → Détecter les fragments → Simuler chaque fragment → Combiner les spectres
```

La condition de coupure : si **|J(A,B) / Δν(A,B)| < ε** (typiquement 0.05), les deux groupes sont indépendants.

---

## Étape 1 : Construire le graphe de couplage

```python
import networkx as nx

def build_coupling_graph(frequencies_hz, J_matrix, threshold_ratio=0.05):
    """
    Construit le graphe de connectivité des spins.
    Une arête existe si |J_ij / delta_nu_ij| > threshold.
    """
    n = len(frequencies_hz)
    G = nx.Graph()
    G.add_nodes_from(range(n))

    for i in range(n):
        for j in range(i + 1, n):
            Jij = abs(J_matrix[i, j])
            if Jij < 0.01:
                continue  # Pas de couplage

            delta_nu = abs(frequencies_hz[i] - frequencies_hz[j])

            # Éviter division par zéro (spins dégénérés = fortement couplés)
            if delta_nu < 0.001:
                G.add_edge(i, j, weight=Jij)
                continue

            ratio = Jij / delta_nu
            if ratio > threshold_ratio:
                G.add_edge(i, j, weight=Jij)

    return G
```

---

## Étape 2 : Détecter les fragments

```python
def find_fragments(frequencies_hz, J_matrix, threshold_ratio=0.05):
    """
    Retourne les composantes connexes = fragments indépendants.

    Returns:
        Liste de listes d'indices, ex: [[0,1,2], [3,4,5,6], [7]]
    """
    G = build_coupling_graph(frequencies_hz, J_matrix, threshold_ratio)

    fragments = list(nx.connected_components(G))
    fragments = [sorted(list(f)) for f in fragments]

    print(f"Système de {len(frequencies_hz)} spins → {len(fragments)} fragments:")
    for i, frag in enumerate(fragments):
        print(f"  Fragment {i+1}: spins {frag} ({len(frag)} spins)")

    return fragments
```

---

## Étape 3 : Simuler chaque fragment séparément

```python
def simulate_fragment(fragment_indices, frequencies_hz, J_matrix,
                      intensity_threshold=0.001):
    """
    Simule un fragment de spins avec votre fonction existante.
    """
    n = len(fragment_indices)

    # Extraire les fréquences du fragment
    freqs_fragment = [frequencies_hz[i] for i in fragment_indices]

    # Extraire la sous-matrice J
    J_fragment = np.zeros((n, n))
    for local_i, global_i in enumerate(fragment_indices):
        for local_j, global_j in enumerate(fragment_indices):
            J_fragment[local_i, local_j] = J_matrix[global_i, global_j]

    if n > 10:
        print(f"  ⚠️  Fragment de {n} spins encore trop grand!")
        print(f"      → Augmenter threshold_ratio ou utiliser Spinach")
        return []

    return simulate_spin_system_exact(freqs_fragment, J_fragment,
                                      intensity_threshold)
```

---

## Étape 4 : Fonction principale de fragmentation

```python
def simulate_large_spin_system(frequencies_hz, J_matrix,
                               threshold_ratio=0.05,
                               intensity_threshold=0.001):
    """
    Simule un grand système de spins par fragmentation.

    Pipeline:
        1. Détecter les fragments faiblement couplés
        2. Simuler chaque fragment exactement
        3. Combiner tous les spectres

    Args:
        frequencies_hz : fréquences de Larmor (Hz)
        J_matrix       : matrice des couplages J (Hz)
        threshold_ratio: seuil |J/Δν| pour couper une arête (défaut 0.05)

    Returns:
        Liste de (freq_hz, intensité) — spectre complet combiné
    """
    n = len(frequencies_hz)
    print(f"\nSimulation de {n} spins par fragmentation")
    print(f"Seuil J/Δν = {threshold_ratio}")
    print("-" * 40)

    # Cas simple: assez petit pour simulation directe
    if n <= 10:
        return simulate_spin_system_exact(frequencies_hz, J_matrix,
                                          intensity_threshold)

    # Détecter les fragments
    fragments = find_fragments(frequencies_hz, J_matrix, threshold_ratio)

    # Simuler chaque fragment
    all_transitions = []

    for i, frag in enumerate(fragments):
        print(f"\nFragment {i+1} ({len(frag)} spins): {frag}")

        transitions = simulate_fragment(frag, frequencies_hz, J_matrix,
                                        intensity_threshold)

        if transitions:
            all_transitions.extend(transitions)
            print(f"  → {len(transitions)} transitions trouvées")

    # Trier par fréquence
    all_transitions.sort(key=lambda x: x[0])

    print(f"\nTotal: {len(all_transitions)} transitions")
    return all_transitions
```

---

## Exemple concret : chaîne octane CH₃-(CH₂)₆-CH₃

```python
def build_alkyl_chain(n_ch2, spectrometer_mhz=300.0):
    """
    Construit les paramètres RMN d'une chaîne -(CH₂)ₙ-CH₃.

    Convention: spins ordonnés le long de la chaîne
        CH₃(a) - CH₂(1) - CH₂(2) - ... - CH₂(n) - CH₃(b)
    """
    # Déplacements chimiques (ppm)
    delta_ch3 = 0.88
    delta_ch2_alpha = 1.26   # CH₂ adjacent au CH₃
    delta_ch2_inner = 1.26   # CH₂ intérieurs (~identiques)

    # Construire la liste des spins (2H par groupe)
    spins = []

    # CH₃ terminal a (3 protons)
    spins += [delta_ch3 * spectrometer_mhz] * 3

    # Groupes CH₂ (2 protons chacun)
    for k in range(n_ch2):
        delta = delta_ch2_alpha if k == 0 else delta_ch2_inner
        spins += [delta * spectrometer_mhz] * 2

    # CH₃ terminal b (3 protons)
    spins += [delta_ch3 * spectrometer_mhz] * 3

    n = len(spins)
    J = np.zeros((n, n))

    # Remplir les couplages (simplifié: par groupe)
    # ³J vicinal entre groupes adjacents ≈ 7 Hz
    # ²J géminal dans CH₂ ≈ -11 Hz

    # ... (logique de couplage selon position)

    return np.array(spins), J


# Test avec octane (8C → 18H)
freqs, J = build_alkyl_chain(n_ch2=6, spectrometer_mhz=300.0)
transitions = simulate_large_spin_system(freqs, J, threshold_ratio=0.05)
```

---

## Résultat attendu pour vos chaînes

| Chaîne    | Spins | Fragments typiques (300 MHz) | Faisable ?                 |
| --------- | ----- | ---------------------------- | -------------------------- |
| -(CH₂)₅-  | 10H   | 1–2 fragments                | ✅ Direct                  |
| -(CH₂)₈-  | 16H   | CH₃ séparés + corps          | ✅ Avec fragmentation      |
| -(CH₂)₁₂- | 24H   | CH₃ × 2 + bloc central 18H   | ⚠️ Bloc central trop grand |

**Le problème restant** : les CH₂ intérieurs ont tous Δδ ≈ 0 → Δν/J ≈ 0 → ils forment **un seul bloc impossible à couper**. C'est la limite fondamentale de la fragmentation — pour ce cas, seul Spinach (RSS) résout vraiment le problème.

Voulez-vous que je complète la fonction `build_alkyl_chain` avec tous les couplages, ou que j'ajoute une correction perturbative pour les couplages négligés ?
