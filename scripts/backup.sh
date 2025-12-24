#!/bin/bash

# ==========================================
# Script de backup Docker NMR
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

BACKUP_DIR="${BACKUP_DIR:-./backups}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
VOLUMES=("nmr_uploads" "nmr_outputs" "nmr_cache")

create_backup() {
    local volume=$1
    local backup_file="$BACKUP_DIR/${volume}_${TIMESTAMP}.tar.gz"
    
    log_info "Sauvegarde de $volume..."
    mkdir -p "$BACKUP_DIR"
    
    docker run --rm \
        -v "$volume:/data:ro" \
        -v "$BACKUP_DIR:/backup" \
        alpine \
        tar czf "/backup/$(basename $backup_file)" -C /data .
    
    local size=$(du -h "$backup_file" | cut -f1)
    log_success "Backup créé: $backup_file ($size)"
}

echo -e "${BLUE}"
cat << "EOF"
╔══════════════════════════════════════════╗
║     Gestion des Backups NMR Docker      ║
╚══════════════════════════════════════════╝
EOF
echo -e "${NC}"

case "${1:-help}" in
    backup)
        if [ -z "$2" ]; then
            log_info "Sauvegarde de tous les volumes..."
            for volume in "${VOLUMES[@]}"; do
                if docker volume ls | grep -q "$volume"; then
                    create_backup "$volume"
                fi
            done
        else
            create_backup "$2"
        fi
        ;;
        
    list)
        log_info "Backups disponibles:"
        echo ""
        if [ -d "$BACKUP_DIR" ]; then
            for volume in "${VOLUMES[@]}"; do
                echo -e "${YELLOW}$volume:${NC}"
                find "$BACKUP_DIR" -name "${volume}_*.tar.gz" | sort -r | head -n 5
                echo ""
            done
        fi
        ;;
        
    help|*)
        echo "Usage: $0 {backup|list|help}"
        echo ""
        echo "Commandes:"
        echo "  backup [volume]   - Sauvegarder un ou tous les volumes"
        echo "  list              - Lister les backups disponibles"
        ;;
esac

echo ""