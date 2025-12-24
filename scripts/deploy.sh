#!/bin/bash

# ==========================================
# Script de déploiement Docker NMR
# ==========================================

set -e

# Couleurs
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Banner
echo -e "${BLUE}"
cat << "EOF"
╔══════════════════════════════════════════╗
║   Déploiement Docker NMR Application    ║
╚══════════════════════════════════════════╝
EOF
echo -e "${NC}"

# Vérifications
log_info "Vérification des prérequis..."

if ! command -v docker &> /dev/null; then
    log_error "Docker n'est pas installé"
    exit 1
fi

if ! command -v docker-compose &> /dev/null && ! docker compose version &> /dev/null; then
    log_error "Docker Compose n'est pas installé"
    exit 1
fi

if docker compose version &> /dev/null; then
    DOCKER_COMPOSE="docker compose"
else
    DOCKER_COMPOSE="docker-compose"
fi

log_success "Docker et Docker Compose sont installes [OK]"

# Verifier .env
if [ ! -f "backend/.env.docker" ]; then
    log_warning "Fichier backend/.env.docker introuvable"
    if [ -f ".env.example" ]; then
        log_info "Copie de .env.example vers backend/.env.docker"
        mkdir -p backend
        cp .env.example backend/.env.docker
        log_warning "[WARN] Editez backend/.env.docker avec vos configurations"
        read -p "Appuyez sur Entree pour continuer..."
    else
        log_error "Aucun fichier .env.example trouvé"
        exit 1
    fi
fi

# Vérifier structure
log_info "Vérification de la structure..."
required_dirs=("backend" "frontend" "docker")
for dir in "${required_dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        log_error "Dossier manquant: $dir"
        exit 1
    fi
done
log_success "Structure valide [OK]"

# Mode
MODE="${1:-build}"

case "$MODE" in
    build)
        log_info "Mode: Build + Start"
        log_info "Construction des images Docker..."
        $DOCKER_COMPOSE build --no-cache
        log_success "Images construites [OK]"
        ;;
    start)
        log_info "Mode: Start"
        ;;
    restart)
        log_info "Mode: Restart"
        $DOCKER_COMPOSE down
        ;;
    stop)
        log_info "Mode: Stop"
        $DOCKER_COMPOSE down
        log_success "Services arrêtés"
        exit 0
        ;;
    clean)
        log_warning "Mode: Clean (suppression totale)"
        read -p "[WARN] Etes-vous sur ? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            $DOCKER_COMPOSE down -v --rmi all
            docker system prune -f
            log_success "Nettoyage complet effectué"
        else
            log_info "Nettoyage annulé"
        fi
        exit 0
        ;;
    *)
        log_error "Mode invalide: $MODE"
        echo "Usage: $0 {build|start|restart|stop|clean}"
        exit 1
        ;;
esac

# Démarrage
log_info "Démarrage des services..."
$DOCKER_COMPOSE up -d

# Attente
log_info "Attente du démarrage..."
sleep 5

# Verification
if docker ps | grep -q nmr-backend; then
    log_success "Backend demarre [OK]"
else
    log_error "Le backend n'a pas demarre"
    $DOCKER_COMPOSE logs backend
    exit 1
fi

if docker ps | grep -q nmr-frontend; then
    log_success "Frontend demarre [OK]"
else
    log_error "Le frontend n'a pas demarre"
    $DOCKER_COMPOSE logs frontend
    exit 1
fi

# Résultat
FRONTEND_PORT=$(grep FRONTEND_PORT .env.example 2>/dev/null | cut -d'=' -f2 || echo "3000")

echo ""
log_success "🎉 Déploiement réussi !"
echo ""
echo -e "${GREEN}╔══════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║         Application démarrée !           ║${NC}"
echo -e "${GREEN}╚══════════════════════════════════════════╝${NC}"
echo ""
echo -e "  📱 Frontend:  ${BLUE}http://localhost:${FRONTEND_PORT:-3000}${NC}"
echo -e "  🔧 Backend:   ${BLUE}http://localhost:4000${NC} (via proxy)"
echo ""
echo -e "${YELLOW}Commandes utiles:${NC}"
echo -e "  • Logs:      ${BLUE}docker-compose logs -f${NC}"
echo -e "  • Arrêter:   ${BLUE}./scripts/deploy.sh stop${NC}"
echo -e "  • Redémarrer:${BLUE}./scripts/deploy.sh restart${NC}"
echo ""