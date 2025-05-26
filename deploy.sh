#!/bin/bash
echo "============================================"
echo "GUIdedRNA Deployment Script"
echo "============================================"

# First, mount Windows drives for WSL2 compatibility
echo "Setting up Windows drive access..."
if [ -f "./mount_drives.sh" ]; then
    chmod +x ./mount_drives.sh
    ./mount_drives.sh
else
    echo "⚠️  mount_drives.sh not found, skipping drive mounting"
    echo "   Windows drives may not be accessible in the container"
fi

echo ""

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed"
    exit 1
fi

# Create necessary directories
echo "Creating directories..."
mkdir -p data output
chmod 755 data output

echo "Building GUIdedRNA Docker image..."

# Build the Docker image
docker build -t guidedrna:latest .

if [ $? -ne 0 ]; then
    echo "❌ Docker build failed!"
    exit 1
fi

echo "✓ Docker build completed successfully"

# Stop and remove existing container if it exists
echo "Cleaning up existing containers..."
docker stop guidedrna-app 2>/dev/null || true
docker rm guidedrna-app 2>/dev/null || true

# Detect available drives
echo "Detecting available Windows drives..."
AVAILABLE_DRIVES=()
VOLUME_MOUNTS=()

# Base volumes that are always included
VOLUME_MOUNTS+=("-v" "$(pwd)/data:/data")
VOLUME_MOUNTS+=("-v" "$(pwd)/output:/output")

# Check for WSL2 mounted drives
if [ -d "/mnt" ]; then
    echo "Found /mnt directory, checking for Windows drives..."
    VOLUME_MOUNTS+=("-v" "/mnt:/host_mnt:ro")
    
    for letter in c d e f g h; do
        if [ -d "/mnt/$letter" ] && [ -r "/mnt/$letter" ]; then
            echo "  ✓ Found accessible drive: $letter"
            AVAILABLE_DRIVES+=("$letter")
            VOLUME_MOUNTS+=("-v" "/mnt/$letter:/mnt/$letter:ro")
        elif [ -d "/mnt/$letter" ]; then
            echo "  ⚠️  Drive $letter exists but may not be accessible"
            # Try to mount it anyway, it might work in Docker
            VOLUME_MOUNTS+=("-v" "/mnt/$letter:/mnt/$letter:ro")
        fi
    done
fi

# Check for home directory
if [ -d "/home" ]; then
    VOLUME_MOUNTS+=("-v" "/home:/host_home:ro")
fi

echo ""
echo "Starting GUIdedRNA container with volume mounts..."
echo "Detected drives: ${AVAILABLE_DRIVES[*]:-none}"

# Run the container with all detected volume mounts
docker run -d \
    -p 3838:3838 \
    "${VOLUME_MOUNTS[@]}" \
    --user "$(id -u):$(id -g)" \
    --name guidedrna-app \
    --restart unless-stopped \
    guidedrna:latest

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ GUIdedRNA is now running!"
    echo "✓ Application URL: http://localhost:3838"
    echo "✓ Data directory: $(pwd)/data"
    echo "✓ Output directory: $(pwd)/output"
    
    if [ ${#AVAILABLE_DRIVES[@]} -gt 0 ]; then
        echo "✓ Windows drives accessible: ${AVAILABLE_DRIVES[*]}"
    fi
    
    echo ""
    echo "To view logs: docker logs guidedrna-app"
    echo "To stop: docker stop guidedrna-app"
    echo "To restart: docker restart guidedrna-app"
else
    echo "❌ Failed to start container!"
    echo "Check the error messages above for details."
    exit 1
fi

echo ""
echo "============================================"
echo "Deployment complete!"
echo "Open your web browser and go to:"
echo "http://localhost:3838"
echo "============================================"