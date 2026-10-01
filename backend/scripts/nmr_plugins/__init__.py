# nmr_plugins/__init__.py
# Modular NMR plugins

from typing import Dict, Any, Optional
import importlib

# Available plugin registry
# Dynamics is temporarily disabled for debugging.
AVAILABLE_PLUGINS = {
    'isotopomers': 'nmr_plugins.isotopomers',
    # 'dynamics': 'nmr_plugins.dynamics',  # Disabled due to J-coupling issues.
    'exotic_carbocycles': 'nmr_plugins.exotic_carbocycles',
    'organometallic': 'nmr_plugins.organometallic',
    'polymers': 'nmr_plugins.polymers',
}

# Loaded module cache
_loaded_plugins = {}

def load_plugin(plugin_name: str):
    """Load a plugin on demand."""
    if plugin_name not in _loaded_plugins:
        if plugin_name in AVAILABLE_PLUGINS:
            module_name = AVAILABLE_PLUGINS[plugin_name]
            _loaded_plugins[plugin_name] = importlib.import_module(module_name)
        else:
            raise ValueError(f"Plugin '{plugin_name}' non disponible")
    return _loaded_plugins[plugin_name]

def get_all_couplings(mol, h1_idx: int, h2_idx: int, **kwargs) -> Dict[str, float]:
    """Try all plugins in priority order.

    Args:
        mol: RDKit molecule.
        h1_idx, h2_idx: Hydrogen atom indices.
        **kwargs: Additional conditions such as temperature and frequency.

    Returns:
        A dictionary mapping plugin names to detected coupling values."""
    results = {}
    
    # Plugin priority order (dynamics disabled)
    priority_order = ['isotopomers', 'exotic_carbocycles', 'organometallic', 'polymers']
    
    for plugin_name in priority_order:
        if plugin_name in AVAILABLE_PLUGINS:
            try:
                plugin = load_plugin(plugin_name)
                if plugin and hasattr(plugin, 'get_coupling'):
                    coupling = plugin.get_coupling(mol, h1_idx, h2_idx, **kwargs)
                    if coupling is not None:
                        # Require a numeric result rather than a complex object.
                        if isinstance(coupling, (int, float)):
                            results[plugin_name] = float(coupling)
                        else:
                            # Try numeric conversion.
                            try:
                                results[plugin_name] = float(coupling)
                            except (ValueError, TypeError):
                                print(f"Impossible de convertir couplage {plugin_name}: {coupling}")
                                continue
            except Exception as e:
                # Log the error and continue with the other plugins.
                print(f"Erreur plugin {plugin_name}: {e}")
                continue
    
    return results


def get_first_coupling(mol, h1_idx: int, h2_idx: int) -> Optional[float]:
    """Return the first detected coupling for compatibility with the legacy API."""
    # Plugin priority order (dynamics disabled)
    priority_order = [
        'isotopomers',      # More specific
        # 'dynamics',       # Disabled
        'exotic_carbocycles',
        'organometallic',
        'polymers',         # More general
    ]
    
    for plugin_name in priority_order:
        try:
            plugin = load_plugin(plugin_name)
            if hasattr(plugin, 'get_coupling'):
                J = plugin.get_coupling(mol, h1_idx, h2_idx)
                if J is not None:
                    return J
        except Exception:
            continue  # Try the next plugin after a failure.
    
    return None

__all__ = ['load_plugin', 'get_all_couplings', 'AVAILABLE_PLUGINS']