#!/bin/bash

echo "==============================================="
echo "    GUIdedRNA - One-Click Installation"
echo "==============================================="
echo
echo "This will install and run GUIdedRNA."
echo "Process may take 5-15 minutes."
echo
echo "Press Enter to continue..."
read

# Create logs directory
mkdir -p logs

# Check if Docker is installed
echo "[1/4] Checking Docker..."
if ! command -v docker &> /dev/null; then
    echo "❌ Docker not found!"
    echo
    echo "Please install Docker first:"
    echo "• Mac: https://docs.docker.com/docker-for-mac/install/"
    echo "• Linux: https://docs.docker.com/engine/install/"
    echo
    echo "After installing Docker:"
    echo "1. Start Docker"
    echo "2. Run this script again"
    echo
    echo "Press Enter to exit..."
    read
    exit 1
fi

# Check if Docker is running
echo "[2/4] Checking if Docker is running..."
if ! docker info &> /dev/null; then
    echo "❌ Docker is installed but not running!"
    echo "Please start Docker and try again."
    echo
    echo "Press Enter to exit..."
    read
    exit 1
fi

echo "✅ Docker is ready!"

# Build the application
echo "[3/4] Building GUIdedRNA (this may take 10-15 minutes)..."
if ! docker build -t guidedrna:latest .; then
    echo "❌ Build failed! Check that Docker is running properly."
    echo "Press Enter to exit..."
    read
    exit 1
fi

echo "✅ Build completed!"

# Start the application
echo "[4/4] Starting GUIdedRNA..."

# Stop any existing container
docker stop guidedrna-app &>/dev/null
docker rm guidedrna-app &>/dev/null

# Find available port
PORT=3838
while netstat -tuln 2>/dev/null | grep -q ":$PORT "; do
    PORT=$((PORT + 1))
    if [ $PORT -gt 3850 ]; then
        echo "❌ No available ports found"
        exit 1
    fi
done

# Start container
docker run -d \
    -p $PORT:3838 \
    -v "$(pwd)/data:/data" \
    -v "$(pwd)/output:/output" \
    --name guidedrna-app \
    guidedrna:latest

echo
echo "✅ GUIdedRNA is now running!"
echo "🌐 Open your browser to: http://localhost:$PORT"
echo
echo "To stop GUIdedRNA: docker stop guidedrna-app"
echo

# Try to open browser
if command -v open &> /dev/null; then
    # macOS
    open "http://localhost:$PORT"
elif command -v xdg-open &> /dev/null; then
    # Linux
    xdg-open "http://localhost:$PORT"
fi

echo "Press Enter to close..."
read
