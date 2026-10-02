#!/bin/bash
echo ""
echo "  Starting NMRScout..."
echo ""

if [ ! -f docker-compose.volumes.yml ]; then
    echo "  First run detected, launching installer..."
    bash "$(dirname "$0")/install.sh"
    exit
fi

echo "  Pulling latest Docker images..."
if ! docker compose -f docker-compose.shareable.yml -f docker-compose.volumes.yml pull; then
    echo "  Failed to pull Docker images. Check Docker Desktop and your internet connection."
    exit 1
fi

docker compose -f docker-compose.shareable.yml -f docker-compose.volumes.yml up -d
if [ $? -ne 0 ]; then
    echo "  Failed to start NMRScout."
    exit 1
fi
echo ""
echo "  NMRScout: http://localhost:3000"
echo ""
