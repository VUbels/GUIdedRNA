#!/bin/bash

echo "============================================"
echo "GUIdedRNA Deployment Script"
echo "============================================"

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed"
    echo "Please install Docker first: https://docs.docker.com/get-docker/"
    exit 1
fi

# Check if docker-compose is installed
if ! command -v docker-compose &> /dev/null; then
    echo "Warning: docker-compose is not installed, using docker build instead"
    USE_COMPOSE=false
else
    USE_COMPOSE=true
fi

# Create necessary directories
echo "Creating directories..."
mkdir -p data output

# Set permissions for directories
chmod 755 data output

echo "Building GUIdedRNA package..."

if [ "$USE_COMPOSE" = true ]; then
    echo "Using docker-compose..."
    docker-compose up --build -d
    
    echo ""
    echo "✓ GUIdedRNA is now running!"
    echo "✓ Application URL: http://localhost:3838"
    echo "✓ Data directory: $(pwd)/data"
    echo "✓ Output directory: $(pwd)/output"
    echo ""
    echo "To view logs: docker-compose logs -f"
    echo "To stop: docker-compose down"
    
else
    echo "Using docker build..."
    docker build -t guidedrna:latest .
    
    echo "Starting container..."
    docker run -d \
        -p 3838:3838 \
        -v "$(pwd)/data:/data" \
        -v "$(pwd)/output:/output" \
        --name guidedrna-app \
        guidedrna:latest
    
    echo ""
    echo "✓ GUIdedRNA is now running!"
    echo "✓ Application URL: http://localhost:3838"
    echo "✓ Data directory: $(pwd)/data"
    echo "✓ Output directory: $(pwd)/output"
    echo ""
    echo "To view logs: docker logs -f guidedrna-app"
    echo "To stop: docker stop guidedrna-app"
    echo "To remove: docker rm guidedrna-app"
fi

echo ""
echo "============================================"
echo "Deployment complete!"
echo "Open your web browser and go to:"
echo "http://localhost:3838"
echo "============================================"
