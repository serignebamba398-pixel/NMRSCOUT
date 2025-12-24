import sys
import os
import importlib.util

# Map command names to script paths (relative to this script)
SCRIPT_MAP = {
    'extract_smiles': 'extract_smiles.py',
    'predict_1H': 'predict_1H_real.py',
    'predict_13C': 'predict_13C_real.py',
    'predict_COSY': 'predict_COSY_real.py',
    'predict_HSQC': 'predict_HSQC_real.py',
    'integrated_prediction': 'integrated_nmr_prediction.py',
    'check_dependencies': 'check_dependencies.py',
    'check_ml_dependencies': 'check_ml_dependencies.py'
}

def run_script(script_name, args):
    if script_name not in SCRIPT_MAP:
        print(f"Error: Unknown script '{script_name}'")
        sys.exit(1)

    script_file = SCRIPT_MAP[script_name]
    
    # In a PyInstaller bundle, we might need to adjust how we find the script
    # But since we are bundling everything, we can import the modules directly
    # if we structure them as modules. 
    # However, the existing scripts are scripts.
    # The easiest way in PyInstaller is to bundle them as data files or 
    # use runpy.
    
    import runpy
    
    # If we are frozen (bundled), the scripts should be in the bundle
    if getattr(sys, 'frozen', False):
        base_path = sys._MEIPASS
        script_path = os.path.join(base_path, 'scripts', script_file)
    else:
        base_path = os.path.dirname(os.path.abspath(__file__))
        script_path = os.path.join(base_path, 'scripts', script_file)

    # Adjust sys.argv so the script thinks it was called directly
    # sys.argv[0] is currently this wrapper script
    # sys.argv[1] is the command name
    # sys.argv[2:] are the arguments for the script
    
    sys.argv = [script_path] + args
    
    # Run the script
    try:
        runpy.run_path(script_path, run_name='__main__')
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        print(f"Error running {script_name}: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: py_entry.py <command> [args...]")
        print("Available commands:", ", ".join(SCRIPT_MAP.keys()))
        sys.exit(1)

    command = sys.argv[1]
    script_args = sys.argv[2:]
    run_script(command, script_args)
