#!/bin/bash

echo "============================================"
echo "WSL2 Drive Mounting Script for GUIdedRNA"
echo "============================================"

# Function to check if a drive is mounted
check_mount() {
    local drive=$1
    if mountpoint -q "/mnt/$drive" 2>/dev/null; then
        echo "✓ Drive $drive: is already mounted at /mnt/$drive"
        return 0
    else
        echo "✗ Drive $drive: not mounted"
        return 1
    fi
}

# Function to mount a drive
mount_drive() {
    local drive=$1
    local mount_point="/mnt/$drive"
    
    echo "Attempting to mount drive $drive..."
    
    # Create mount point if it doesn't exist
    sudo mkdir -p "$mount_point"
    
    # Try to mount the drive with appropriate permissions
    if sudo mount -t drvfs "${drive}:" "$mount_point" -o uid=$(id -u),gid=$(id -g),metadata,umask=022,fmask=111; then
        echo "✓ Successfully mounted drive $drive at $mount_point"
        return 0
    else
        echo "✗ Failed to mount drive $drive"
        return 1
    fi
}

# Function to unmount a drive (for cleanup if needed)
unmount_drive() {
    local drive=$1
    local mount_point="/mnt/$drive"
    
    if mountpoint -q "$mount_point" 2>/dev/null; then
        echo "Unmounting drive $drive..."
        sudo umount "$mount_point"
        echo "✓ Drive $drive unmounted"
    fi
}

# Check current WSL configuration
echo ""
echo "Current WSL mount configuration:"
if [ -f /etc/wsl.conf ]; then
    echo "Found /etc/wsl.conf:"
    cat /etc/wsl.conf
else
    echo "No /etc/wsl.conf found. Creating optimal configuration..."
    
    # Create optimal WSL configuration
    sudo tee /etc/wsl.conf > /dev/null << 'EOF'
[automount]
enabled = true
root = /mnt/
options = "metadata,uid=1000,gid=1000,umask=022,fmask=111"
mountFsTab = true

[interop]
enabled = true
appendWindowsPath = true

[user]
default = $(whoami)
EOF
    
    echo "✓ Created /etc/wsl.conf with optimal settings"
    echo "⚠️  You may need to restart WSL for these changes to take effect:"
    echo "   Run 'wsl --shutdown' from Windows PowerShell, then restart WSL"
fi

echo ""
echo "Checking and mounting Windows drives..."

# List of common drive letters to check
DRIVES=("c" "d" "e" "f" "g" "h")

# Track which drives we successfully mount
MOUNTED_DRIVES=()

for drive in "${DRIVES[@]}"; do
    echo ""
    echo "Checking drive $drive..."
    
    # Check if the Windows drive exists (try to access it from Windows)
    if [ -e "/mnt/$drive" ] || [ -L "/mnt/$drive" ]; then
        if check_mount "$drive"; then
            MOUNTED_DRIVES+=("$drive")
        else
            if mount_drive "$drive"; then
                MOUNTED_DRIVES+=("$drive")
            fi
        fi
    else
        echo "✗ Drive $drive: not found in Windows"
    fi
done

echo ""
echo "============================================"
echo "Drive Mounting Summary"
echo "============================================"

if [ ${#MOUNTED_DRIVES[@]} -eq 0 ]; then
    echo "⚠️  No drives were successfully mounted"
    echo ""
    echo "Troubleshooting suggestions:"
    echo "1. Restart WSL completely:"
    echo "   - Run 'wsl --shutdown' from Windows PowerShell"
    echo "   - Restart your WSL Ubuntu session"
    echo "   - Run this script again"
    echo ""
    echo "2. Check Windows drives manually:"
    echo "   - Open Windows File Explorer"
    echo "   - Note which drive letters exist (C:, D:, etc.)"
    echo "   - Ensure the drives are accessible"
    echo ""
    echo "3. Try manual mounting:"
    echo "   sudo mkdir -p /mnt/c"
    echo "   sudo mount -t drvfs C: /mnt/c"
else
    echo "✓ Successfully mounted drives: ${MOUNTED_DRIVES[*]}"
    echo ""
    echo "Mounted drive locations:"
    for drive in "${MOUNTED_DRIVES[@]}"; do
        echo "  Windows $drive: -> /mnt/$drive"
        # Show some content to verify
        echo "    Sample content: $(ls /mnt/$drive 2>/dev/null | head -3 | tr '\n' ' ')..."
    done
fi

echo ""
echo "Current mount points:"
df -h | grep "/mnt" || echo "No /mnt drives currently mounted"

echo ""
echo "============================================"
echo "Mount script completed!"
echo "============================================"
