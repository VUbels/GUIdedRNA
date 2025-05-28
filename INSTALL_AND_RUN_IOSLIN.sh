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

# Create data and output directories if they don't exist
mkdir -p data output

# Detect operating system and set up appropriate volume mounts
OS_TYPE=$(uname -s)
VOLUME_MOUNTS=""

# Common mounts for all Unix systems
VOLUME_MOUNTS+="-v $(pwd)/data:/data "
VOLUME_MOUNTS+="-v $(pwd)/output:/output "
VOLUME_MOUNTS+="-v $HOME:/host_home "

case "$OS_TYPE" in
    "Linux")
        echo "🐧 Detected Linux - mounting system drives..."
        VOLUME_MOUNTS+="-v /:/host_system:ro "
        VOLUME_MOUNTS+="-v /media:/host_media:ro "
        VOLUME_MOUNTS+="-v /mnt:/host_mnt:ro "
        
        # Check for WSL and Windows drives
        if grep -qi microsoft /proc/version 2>/dev/null; then
            echo "🪟 WSL detected - mounting Windows drives..."
            for drive in c d e f g h i j; do
                if [ -d "/mnt/$drive" ]; then
                    VOLUME_MOUNTS+="-v /mnt/$drive:/host_drives/${drive^^}:ro "
                fi
            done
        fi
        ;;
    "Darwin")
        echo "🍎 Detected macOS - mounting system drives..."
        VOLUME_MOUNTS+="-v /:/host_system:ro "
        VOLUME_MOUNTS+="-v /Volumes:/host_volumes:ro "
        ;;
    *)
        echo "❓ Unknown OS - using basic mounts..."
        VOLUME_MOUNTS+="-v /:/host_system:ro "
        ;;
esac

# Start container with comprehensive volume mounts
echo "Starting container with volume mounts..."
docker run -d \
    -p $PORT:3838 \
    $VOLUME_MOUNTS \
    --name guidedrna-app \
    guidedrna:latest

echo
echo "✅ GUIdedRNA is now running!"
echo "🌐 Open your browser to: http://localhost:$PORT"
echo
echo "📁 Available drives should now be visible in the file browser:"
case "$OS_TYPE" in
    "Linux")
        echo "🏠 Your home folder: 'Host Home'"
        echo "💾 System drives: 'Host System', 'Host Media', 'Host Mnt'"
        if grep -qi microsoft /proc/version 2>/dev/null; then
            echo "🪟 Windows drives: 'C:', 'D:', 'E:', etc."
        fi
        ;;
    "Darwin")
        echo "🏠 Your home folder: 'Host Home'"
        echo "💾 External drives: 'Host Volumes'"
        echo "🖥️  System: 'Host System'"
        ;;
esac
echo
echo "To stop GUIdedRNA: docker stop guidedrna-app"
echo

# Wait for container to start
echo "Waiting for application to start..."
sleep 10

# Try to open browser
if command -v open &> /dev/null; then
    # macOS
    open "http://localhost:$PORT"
elif command -v xdg-open &> /dev/null; then
    # Linux
    xdg-open "http://localhost:$PORT"
fi

echo
echo "If the browser doesn't open automatically, manually go to:"
echo "http://localhost:$PORT"
echo
echo "Press Enter to close..."
read