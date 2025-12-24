# verify_dependencies.py - Vérification complète des dépendances du système NMR ML
import sys
import os
import subprocess
import importlib
import platform
from pathlib import Path

class DependencyChecker:
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.success_count = 0
        self.total_checks = 0
        
    def log_success(self, message):
        print(f"✅ {message}")
        self.success_count += 1
    
    def log_warning(self, message):
        print(f"⚠️  {message}")
        self.warnings.append(message)
    
    def log_error(self, message):
        print(f"❌ {message}")
        self.errors.append(message)
    
    def check_python_version(self):
        """Vérifier la version de Python"""
        print("\n🔍 Vérification de Python...")
        self.total_checks += 1
        
        version = sys.version_info
        version_str = f"{version.major}.{version.minor}.{version.micro}"
        
        if version.major == 3 and version.minor >= 8:
            self.log_success(f"Python {version_str} (compatible)")
        elif version.major == 3 and version.minor >= 6:
            self.log_warning(f"Python {version_str} (fonctionnel mais Python 3.8+ recommandé)")
        else:
            self.log_error(f"Python {version_str} (version 3.8+ requise)")
    
    def check_python_packages(self):
        """Vérifier les packages Python critiques"""
        print("\n🔍 Vérification des packages Python...")
        
        critical_packages = [
            ('torch', 'PyTorch pour ML'),
            ('rdkit', 'RDKit pour chimie'),
            ('pandas', 'Manipulation de données'),
            ('numpy', 'Calculs numériques'),
        ]
        
        optional_packages = [
            ('dgl', 'Deep Graph Library'),
            ('dgllife', 'DGL pour sciences de la vie'),
            ('sklearn', 'Scikit-learn'),
            ('matplotlib', 'Visualisation'),
            ('tqdm', 'Barres de progression'),
        ]
        
        # Packages critiques
        for package, description in critical_packages:
            self.total_checks += 1
            try:
                if package == 'rdkit':
                    # RDKit a un nom d'import différent
                    importlib.import_module('rdkit.Chem')
                else:
                    importlib.import_module(package)
                
                # Obtenir la version si possible
                try:
                    if package == 'torch':
                        import torch
                        version = torch.__version__
                    elif package == 'rdkit':
                        import rdkit
                        version = rdkit.__version__
                    elif package == 'pandas':
                        import pandas
                        version = pandas.__version__
                    elif package == 'numpy':
                        import numpy
                        version = numpy.__version__
                    else:
                        version = "installé"
                    
                    self.log_success(f"{package} {version} - {description}")
                except:
                    self.log_success(f"{package} - {description}")
                    
            except ImportError:
                self.log_error(f"{package} manquant - {description}")
        
        # Packages optionnels
        for package, description in optional_packages:
            self.total_checks += 1
            try:
                importlib.import_module(package)
                self.log_success(f"{package} - {description}")
            except ImportError:
                self.log_warning(f"{package} manquant (optionnel) - {description}")
    
    def check_torch_functionality(self):
        """Vérifier que PyTorch fonctionne correctement"""
        print("\n🔍 Test de PyTorch...")
        self.total_checks += 1
        
        try:
            import torch
            
            # Test création de tensor
            x = torch.randn(2, 3)
            y = torch.randn(3, 2)
            z = torch.mm(x, y)
            
            # Vérifier CUDA (optionnel)
            cuda_available = torch.cuda.is_available()
            if cuda_available:
                self.log_success("PyTorch avec support CUDA")
            else:
                self.log_success("PyTorch (CPU uniquement)")
            
            self.log_success("Test PyTorch réussi")
            
        except Exception as e:
            self.log_error(f"Erreur PyTorch: {str(e)}")
    
    def check_rdkit_functionality(self):
        """Vérifier que RDKit fonctionne correctement"""
        print("\n🔍 Test de RDKit...")
        self.total_checks += 1
        
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Test création molécule
            mol = Chem.MolFromSmiles('CCO')
            if mol is not None:
                # Test calcul propriétés
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                
                self.log_success("Test RDKit réussi")
            else:
                self.log_error("RDKit ne peut pas créer de molécules")
                
        except Exception as e:
            self.log_error(f"Erreur RDKit: {str(e)}")
    
    def check_node_npm(self):
        """Vérifier Node.js et npm"""
        print("\n🔍 Vérification de Node.js et npm...")
        
        # Node.js
        self.total_checks += 1
        try:
            result = subprocess.run(['node', '--version'], 
                                  capture_output=True, text=True, check=True)
            node_version = result.stdout.strip()
            
            # Extraire le numéro de version
            version_num = int(node_version.replace('v', '').split('.')[0])
            if version_num >= 16:
                self.log_success(f"Node.js {node_version} (compatible)")
            else:
                self.log_warning(f"Node.js {node_version} (version 16+ recommandée)")
                
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.log_error("Node.js non trouvé")
        
        # npm
        self.total_checks += 1
        try:
            result = subprocess.run(['npm', '--version'], 
                                  capture_output=True, text=True, check=True)
            npm_version = result.stdout.strip()
            self.log_success(f"npm {npm_version}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.log_error("npm non trouvé")
    
    def check_project_structure(self):
        """Vérifier la structure du projet"""
        print("\n🔍 Vérification de la structure du projet...")
        
        required_dirs = [
            'backend',
            'backend/models',
            'backend/scripts', 
            'backend/uploads',
            'backend/outputs',
            'frontend',
        ]
        
        recommended_dirs = [
            'backend/ml_dependencies',
            'backend/cache',
            'installation',
        ]
        
        # Dossiers requis
        for dir_path in required_dirs:
            self.total_checks += 1
            if os.path.exists(dir_path):
                self.log_success(f"Dossier {dir_path}")
            else:
                self.log_error(f"Dossier manquant: {dir_path}")
        
        # Dossiers recommandés
        for dir_path in recommended_dirs:
            self.total_checks += 1
            if os.path.exists(dir_path):
                self.log_success(f"Dossier {dir_path}")
            else:
                self.log_warning(f"Dossier recommandé manquant: {dir_path}")
    
    def check_models(self):
        """Vérifier la présence des modèles ML"""
        print("\n🔍 Vérification des modèles ML...")
        
        models_dir = Path('backend/models')
        if not models_dir.exists():
            self.log_error("Dossier backend/models manquant")
            return
        
        required_models = [
            '1H_sparsified_proposed_proposed_set2set_1.pt',
            '13C_sparsified_proposed_proposed_set2set_1.pt'
        ]
        
        for model_file in required_models:
            self.total_checks += 1
            model_path = models_dir / model_file
            
            if model_path.exists():
                size_mb = model_path.stat().st_size / (1024 * 1024)
                self.log_success(f"Modèle {model_file} ({size_mb:.1f} MB)")
            else:
                self.log_error(f"Modèle manquant: {model_file}")
    
    def check_ml_dependencies(self):
        """Vérifier les dépendances ML personnalisées"""
        print("\n🔍 Vérification des dépendances ML...")
        
        ml_dir = Path('backend/ml_dependencies')
        if not ml_dir.exists():
            self.log_warning("Dossier ml_dependencies manquant")
            return
        
        required_files = [
            'gnn_models/__init__.py',
            'gnn_models/mpnn_proposed.py',
            'model.py',
            'dataset.py',
            'util.py'
        ]
        
        for file_path in required_files:
            self.total_checks += 1
            full_path = ml_dir / file_path
            
            if full_path.exists():
                self.log_success(f"Fichier ML {file_path}")
            else:
                self.log_warning(f"Fichier ML manquant: {file_path}")
    
    def check_configuration(self):
        """Vérifier les fichiers de configuration"""
        print("\n🔍 Vérification de la configuration...")
        
        config_files = [
            ('.env', 'Variables d\'environnement'),
            ('requirements.txt', 'Dépendances Python'),
            ('backend/package.json', 'Configuration backend'),
            ('frontend/package.json', 'Configuration frontend'),
        ]
        
        for file_path, description in config_files:
            self.total_checks += 1
            if os.path.exists(file_path):
                self.log_success(f"{description}: {file_path}")
            else:
                self.log_warning(f"{description} manquant: {file_path}")
    
    def check_system_resources(self):
        """Vérifier les ressources système"""
        print("\n🔍 Vérification des ressources système...")
        
        # Mémoire disponible
        try:
            import psutil
            memory = psutil.virtual_memory()
            memory_gb = memory.total / (1024**3)
            
            if memory_gb >= 8:
                self.log_success(f"Mémoire: {memory_gb:.1f} GB (suffisant)")
            elif memory_gb >= 4:
                self.log_warning(f"Mémoire: {memory_gb:.1f} GB (minimal)")
            else:
                self.log_warning(f"Mémoire: {memory_gb:.1f} GB (peut être insuffisant)")
                
        except ImportError:
            self.log_warning("psutil non disponible, impossible de vérifier la mémoire")
        
        # Espace disque
        try:
            disk_usage = os.statvfs('.')
            free_space_gb = (disk_usage.f_frsize * disk_usage.f_bavail) / (1024**3)
            
            if free_space_gb >= 5:
                self.log_success(f"Espace disque libre: {free_space_gb:.1f} GB")
            else:
                self.log_warning(f"Espace disque libre: {free_space_gb:.1f} GB (faible)")
                
        except (AttributeError, OSError):
            # Windows ou erreur
            self.log_warning("Impossible de vérifier l'espace disque")
    
    def run_all_checks(self):
        """Exécuter tous les tests"""
        print("🔍 VÉRIFICATION DES DÉPENDANCES SYSTÈME NMR ML")
        print("=" * 50)
        
        print(f"Plateforme: {platform.system()} {platform.release()}")
        print(f"Architecture: {platform.machine()}")
        print(f"Python: {sys.version}")
        
        # Exécuter tous les tests
        self.check_python_version()
        self.check_python_packages()
        self.check_torch_functionality()
        self.check_rdkit_functionality()
        self.check_node_npm()
        self.check_project_structure()
        self.check_models()
        self.check_ml_dependencies()
        self.check_configuration()
        self.check_system_resources()
        
        # Résumé final
        print("\n" + "=" * 50)
        print("📊 RÉSUMÉ DE LA VÉRIFICATION")
        print("=" * 50)
        
        success_rate = (self.success_count / self.total_checks) * 100 if self.total_checks > 0 else 0
        
        print(f"✅ Tests réussis: {self.success_count}/{self.total_checks} ({success_rate:.1f}%)")
        print(f"⚠️  Avertissements: {len(self.warnings)}")
        print(f"❌ Erreurs: {len(self.errors)}")
        
        if self.errors:
            print("\n🚨 ERREURS CRITIQUES:")
            for error in self.errors:
                print(f"  • {error}")
            print("\n💡 Actions requises:")
            print("  1. Corrigez les erreurs ci-dessus")
            print("  2. Relancez la vérification")
        
        if self.warnings:
            print("\n⚠️  AVERTISSEMENTS:")
            for warning in self.warnings:
                print(f"  • {warning}")
            print("\n💡 Recommandations:")
            print("  1. Installez les composants manquants pour une expérience optimale")
            print("  2. Consultez la documentation d'installation")
        
        if not self.errors and success_rate >= 80:
            print("\n🎉 SYSTÈME PRÊT!")
            print("Le système NMR ML devrait fonctionner correctement.")
            print("\nPour démarrer:")
            print("  • Backend: cd backend && npm start")
            print("  • Frontend: cd frontend && npm start")
            print("  • Ou utilisez: ./start_system.sh (Linux/Mac) ou start_system.bat (Windows)")
            return True
        else:
            print("\n⚠️  SYSTÈME NON PRÊT")
            print("Corrigez les erreurs avant de continuer.")
            return False

def main():
    """Fonction principale"""
    checker = DependencyChecker()
    success = checker.run_all_checks()
    
    # Code de sortie pour scripts automatisés
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()