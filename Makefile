# ==========================================
# Makefile pour NMR Docker Application
# ==========================================

.PHONY: help build start stop restart logs clean dev

# Couleurs
BLUE := \033[0;34m
GREEN := \033[0;32m
YELLOW := \033[1;33m
NC := \033[0m

COMPOSE := docker-compose
PROJECT_NAME := nmr-app

help:
	@echo "$(BLUE)╔══════════════════════════════════════════╗$(NC)"
	@echo "$(BLUE)║   NMR Docker Application - Makefile     ║$(NC)"
	@echo "$(BLUE)╚══════════════════════════════════════════╝$(NC)"
	@echo ""
	@echo "$(GREEN)Production:$(NC)"
	@echo "  make setup         - Configuration initiale"
	@echo "  make build         - Build les images"
	@echo "  make start         - Démarrer"
	@echo "  make stop          - Arrêter"
	@echo "  make restart       - Redémarrer"
	@echo "  make logs          - Voir les logs"
	@echo "  make status        - État des services"
	@echo "  make clean         - Nettoyer tout"
	@echo ""

setup:
	@echo "$(BLUE)[SETUP]$(NC) Configuration initiale..."
	@if [ ! -f backend/.env.docker ]; then \
		mkdir -p backend; \
		cp .env.example backend/.env.docker; \
		echo "$(GREEN)[OK]$(NC) Fichier .env.docker créé"; \
	fi
	@chmod +x scripts/*.sh
	@echo "$(GREEN)[SUCCESS]$(NC) Configuration terminée !"

build: setup
	@echo "$(BLUE)[BUILD]$(NC) Construction des images..."
	$(COMPOSE) build --no-cache
	@echo "$(GREEN)[SUCCESS]$(NC) Images construites !"

start:
	@echo "$(BLUE)[START]$(NC) Démarrage..."
	$(COMPOSE) up -d
	@sleep 3
	@$(COMPOSE) ps
	@echo ""
	@echo "$(GREEN)Application démarrée !$(NC)"
	@echo "  Frontend: http://localhost:3000"

stop:
	@echo "$(BLUE)[STOP]$(NC) Arrêt..."
	$(COMPOSE) down

restart: stop start

logs:
	$(COMPOSE) logs -f

status:
	@$(COMPOSE) ps

clean:
	@echo "$(YELLOW)[WARNING]$(NC) Suppression totale !"
	@read -p "Êtes-vous sûr ? (y/N): " -n 1 -r; \
	echo; \
	if [[ $$REPLY =~ ^[Yy]$$ ]]; then \
		$(COMPOSE) down -v --rmi all; \
		docker system prune -f; \
		echo "$(GREEN)[SUCCESS]$(NC) Nettoyage terminé"; \
	fi