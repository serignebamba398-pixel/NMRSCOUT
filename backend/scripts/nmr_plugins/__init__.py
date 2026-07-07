# nmr_plugins/__init__.py
# Système modulaire de plugins NMR pour couverture complète

from typing import Dict, Any, Optional
import importlib

# Registry des plugins disponibles
# NOTE: dynamics temporairement désactivé pour debug
AVAILABLE_PLUGINS = {
    'isotopomers': 'nmr_plugins.isotopomers',
    # 'dynamics': 'nmr_plugins.dynamics',  # DÉSACTIVÉ - cause des problèmes de J coupling
    'exotic_carbocycles': 'nmr_plugins.exotic_carbocycles',
    'organometallic': 'nmr_plugins.organometallic',
    'polymers': 'nmr_plugins.polymers',
}

# Cache des modules chargés
_loaded_plugins = {}

def load_plugin(plugin_name: str):
    """Charge un plugin de façon lazy."""
    if plugin_name not in _loaded_plugins:
        if plugin_name in AVAILABLE_PLUGINS:
            module_name = AVAILABLE_PLUGINS[plugin_name]
            _loaded_plugins[plugin_name] = importlib.import_module(module_name)
        else:
            raise ValueError(f"Plugin '{plugin_name}' non disponible")
    return _loaded_plugins[plugin_name]

def get_all_couplings(mol, h1_idx: int, h2_idx: int, **kwargs) -> Dict[str, float]:
    """
    Orchestrateur principal - teste tous les plugins par ordre de priorité.
    
    Args:
        mol: Molécule RDKit
        h1_idx, h2_idx: Indices des atomes d'hydrogène
        **kwargs: Arguments additionnels (température, fréquence, etc.)
        
    Returns:
        Dict {plugin_name: coupling_value} des couplages détectés
    """
    results = {}
    
    # Ordre de priorité des plugins (dynamics désactivé)
    priority_order = ['isotopomers', 'exotic_carbocycles', 'organometallic', 'polymers']
    
    for plugin_name in priority_order:
        if plugin_name in AVAILABLE_PLUGINS:
            try:
                plugin = load_plugin(plugin_name)
                if plugin and hasattr(plugin, 'get_coupling'):
                    coupling = plugin.get_coupling(mol, h1_idx, h2_idx, **kwargs)
                    if coupling is not None:
                        # S'assurer que c'est un nombre, pas un objet complexe
                        if isinstance(coupling, (int, float)):
                            results[plugin_name] = float(coupling)
                        else:
                            # Essayer de convertir
                            try:
                                results[plugin_name] = float(coupling)
                            except (ValueError, TypeError):
                                print(f"Impossible de convertir couplage {plugin_name}: {coupling}")
                                continue
            except Exception as e:
                # Log l'erreur mais continue avec les autres plugins
                print(f"Erreur plugin {plugin_name}: {e}")
                continue
    
    return results


def get_first_coupling(mol, h1_idx: int, h2_idx: int) -> Optional[float]:
    """
    Version simple - retourne le premier couplage trouvé.
    Fonction de compatibilité avec l'ancienne API.
    """
    # Ordre de priorité des plugins (dynamics désactivé)
    priority_order = [
        'isotopomers',      # Plus spécifique
        # 'dynamics',       # DÉSACTIVÉ
        'exotic_carbocycles',
        'organometallic',
        'polymers',         # Plus général
    ]
    
    for plugin_name in priority_order:
        try:
            plugin = load_plugin(plugin_name)
            if hasattr(plugin, 'get_coupling'):
                J = plugin.get_coupling(mol, h1_idx, h2_idx)
                if J is not None:
                    return J
        except Exception:
            continue  # Plugin échoue, essayer le suivant
    
    return None

__all__ = ['load_plugin', 'get_all_couplings', 'AVAILABLE_PLUGINS']