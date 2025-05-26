#!/bin/bash

echo "============================================"
echo "GUIdedRNA Deployment Script"
echo "============================================"

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed"
    exit 1
fi

# Create necessary directories
echo "Creating directories..."
mkdir -p data output
chmod 755 data output

echo "Building GUIdedRNA from GitHub..."

# Build and run with docker-compose
if command -v docker-compose &> /dev/null; then
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
    
    docker run -d \
        -p 3838:3838 \
        -v "$(pwd)/data:/data" \
        -v "$(pwd)/output:/output" \
        --name guidedrna-app \
        guidedrna:latest
    
    echo ""
    echo "✓ GUIdedRNA is now running!"
    echo "✓ Application URL: http://localhost:3838"
fi

echo ""
echo "============================================"
echo "Deployment complete!"
echo "Open your web browser and go to:"
echo "http://localhost:3838"
echo "============================================"