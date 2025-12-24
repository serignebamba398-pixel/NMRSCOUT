#!/bin/bash

# ==========================================
# Script de monitoring Docker NMR
# ==========================================

# Couleurs
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

if docker compose version &> /dev/null 2>&1; then
    DOCKER_COMPOSE="docker compose"
else
    DOCKER_COMPOSE="docker-compose"
fi

print_header() {
    echo -e "${BLUE}╔══════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║   Monitoring NMR Docker Application     ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════╝${NC}"
    echo ""
}

check_service_status() {
    local service=$1
    if docker ps --format '{{.Names}}' | grep -q "$service"; then
        echo -e "${GREEN}[OK]${NC} $service est en cours d'execution"
        return 0
    else
        echo -e "${RED}[X]${NC} $service n'est pas en cours d'execution"
        return 1
    fi
}

get_container_stats() {
    local service=$1
    echo -e "\n${YELLOW}Statistiques de $service:${NC}"
    docker stats --no-stream --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}\t{{.NetIO}}" | grep "$service"
}

print_header

case "${1:-status}" in
    status)
        echo -e "${BLUE}État des services:${NC}\n"
        check_service_status "nmr-backend"
        echo ""
        check_service_status "nmr-frontend"
        ;;
        
    stats)
        echo -e "${BLUE}Statistiques:${NC}"
        get_container_stats "nmr-backend"
        get_container_stats "nmr-frontend"
        ;;
        
    logs)
        service="${2:-backend}"
        lines="${3:-50}"
        echo -e "${YELLOW}Logs $service (dernières $lines lignes):${NC}"
        $DOCKER_COMPOSE logs --tail=$lines $service
        ;;
        
    watch)
        echo -e "${YELLOW}Surveillance (Ctrl+C pour arrêter)${NC}\n"
        watch -n 5 "docker stats --no-stream nmr-backend nmr-frontend"
        ;;
        
    help|*)
        echo "Usage: $0 {status|stats|logs|watch|help}"
        echo ""
        echo "Commandes:"
        echo "  status   - État des services"
        echo "  stats    - Statistiques CPU/RAM"
        echo "  logs     - Logs (ex: logs backend 100)"
        echo "  watch    - Surveillance en temps réel"
        ;;
esac

echo ""