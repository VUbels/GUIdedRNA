#!/bin/bash

echo "============================================"
echo "GUIdedRNA Setup Script"
echo "============================================"
echo ""

# Check if running in WSL
if grep -qi microsoft /proc/version 2>/dev/null; then
    echo "✓ Detected WSL2 environment"
    WSL_DETECTED=true
else
    echo "ℹ️  Not running in WSL2"
    WSL_DETECTED=false
fi

# Check Docker installation
if command -v docker &> /dev/null; then
    echo "✓ Docker is installed"
    
    # Check if Docker daemon is running
    if docker ps &> /dev/null; then
        echo "✓ Docker daemon is running"
    else
        echo "❌ Docker daemon is not running or accessible"
        echo "   Please start Docker or check permissions"
        exit 1
    fi
else
    echo "❌ Docker is not installed"
    echo "   Please install Docker first: https://docs.docker.com/get-docker/"
    exit 1
fi

# Check user permissions for Docker
if groups $USER | grep -q docker; then
    echo "✓ User has Docker permissions"
else
    echo "⚠️  User may not have Docker permissions"
    echo "   If you get permission errors, run:"
    echo "   sudo usermod -aG docker $USER"
    echo "   Then log out and log back in"
fi

# Check WSL drive access if in WSL
if [ "$WSL_DETECTED" = true ]; then
    echo ""
    echo "Checking WSL drive access..."
    
    if [ -d "/mnt/c" ]; then
        if [ -r "/mnt/c" ]; then
            echo "✓ C: drive is accessible"
        else
            echo "⚠️  C: drive exists but may not be readable"
        fi
    else
        echo "❌ C: drive not found"
        echo "   WSL may need to be restarted: wsl --shutdown"
    fi
    
    if [ -d "/mnt/d" ]; then
        if [ -r "/mnt/d" ]; then
            echo "✓ D: drive is accessible"
        else
            echo "⚠️  D: drive exists but may not be readable"
        fi
    else
        echo "ℹ️  D: drive not found (this is normal if you don't have a D: drive)"
    fi
fi

echo ""
echo "============================================"
echo "Setup Check Complete!"
echo "============================================"
echo ""

if [ "$WSL_DETECTED" = true ]; then
    echo "For WSL2 users:"
    echo "1. Make sure all your drives are accessible in /mnt/"
    echo "2. If you have issues, try restarting WSL:"
    echo "   - From Windows PowerShell: wsl --shutdown"
    echo "   - Then restart your WSL session"
    echo ""
fi

echo "To deploy GUIdedRNA, run:"
echo "  ./deploy.sh"
echo ""
echo "The deployment script will automatically:"
echo "- Build the Docker image"
echo "- Detect and mount available drives"
echo "- Start the application on http://localhost:3838"
