#!/bin/bash
echo ""
echo "  Starting NMRScout..."
echo ""

if [ ! -f docker-compose.volumes.yml ]; then
    echo "  First run detected, launching installer..."
    bash "$(dirname "$0")/install.sh"
    exit
fi

docker compose -f docker-compose.shareable.yml -f docker-compose.volumes.yml up -d
echo ""
echo "  NMRScout: http://localhost:3000"
echo ""
