#!/bin/bash
# Import a competitor submission file (zip/tar/directory) and organize it
#
# Usage:
#   tools/import_submission.sh submission_file.zip
#   tools/import_submission.sh submission_001_username.tar.gz
#   tools/import_submission.sh submission_directory/
#
# Extracts problem number from filename (e.g., "001" -> problem-001-dimer-ksat)
# Extracts username from filename or submission.json
# Organizes into: submissions/problem-XXX/username/

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"

SUBMISSION_INPUT="${1}"

if [ -z "$SUBMISSION_INPUT" ]; then
    echo "Usage: tools/import_submission.sh <submission_file_or_dir>"
    echo ""
    echo "Examples:"
    echo "  tools/import_submission.sh submission_001_alice.zip"
    echo "  tools/import_submission.sh submission_001_alice.tar.gz"
    echo "  tools/import_submission.sh submission_directory/"
    exit 1
fi

# Check if input exists
if [ ! -e "$SUBMISSION_INPUT" ]; then
    echo "Error: Submission file/directory not found: $SUBMISSION_INPUT"
    exit 1
fi

echo "=========================================="
echo "Importing Submission"
echo "=========================================="
echo "Input: $SUBMISSION_INPUT"
echo ""

# Extract problem number from filename (look for 001, 002, etc.)
PROBLEM_NUM=$(echo "$(basename "$SUBMISSION_INPUT")" | grep -oE '[0-9]{3}' | head -1)

if [ -z "$PROBLEM_NUM" ]; then
    echo "Error: Could not extract problem number from filename."
    echo "Filename should contain problem number like '001', '002', etc."
    echo "Example: submission_001_username.zip"
    exit 1
fi

# Map problem number to problem ID
case "$PROBLEM_NUM" in
    "001")
        PROBLEM_ID="problem-001-dimer-ksat"
        ;;
    "002")
        PROBLEM_ID="problem-002-octahedron"
        ;;
    "003")
        PROBLEM_ID="problem-003-hamiltonian-path"
        ;;
    *)
        echo "Warning: Unknown problem number $PROBLEM_NUM, using problem-$PROBLEM_NUM"
        PROBLEM_ID="problem-$PROBLEM_NUM"
        ;;
esac

echo "Detected problem: $PROBLEM_ID"

# Create temp extraction directory
TEMP_DIR=$(mktemp -d -t submission_import_XXXXXX)
trap "rm -rf $TEMP_DIR" EXIT

# Extract or copy submission
if [ -d "$SUBMISSION_INPUT" ]; then
    echo "Copying directory..."
    cp -r "$SUBMISSION_INPUT"/* "$TEMP_DIR/"
elif [[ "$SUBMISSION_INPUT" == *.zip ]]; then
    echo "Extracting ZIP file..."
    unzip -q "$SUBMISSION_INPUT" -d "$TEMP_DIR"
elif [[ "$SUBMISSION_INPUT" == *.tar.gz ]] || [[ "$SUBMISSION_INPUT" == *.tgz ]]; then
    echo "Extracting TAR.GZ file..."
    tar -xzf "$SUBMISSION_INPUT" -C "$TEMP_DIR"
elif [[ "$SUBMISSION_INPUT" == *.tar ]]; then
    echo "Extracting TAR file..."
    tar -xf "$SUBMISSION_INPUT" -C "$TEMP_DIR"
else
    echo "Error: Unsupported file format. Use .zip, .tar.gz, .tar, or directory"
    exit 1
fi

# Try to extract username from submission.json or filename
USERNAME=""
if [ -f "$TEMP_DIR/submission.json" ]; then
    USERNAME=$(python3 -c "import json; f=open('$TEMP_DIR/submission.json'); d=json.load(f); print(d.get('team_name', d.get('username', '')))" 2>/dev/null || echo "")
fi

# Fallback: extract from filename
if [ -z "$USERNAME" ]; then
    BASENAME=$(basename "$SUBMISSION_INPUT" .zip)
    BASENAME=$(basename "$BASENAME" .tar.gz)
    BASENAME=$(basename "$BASENAME" .tar)
    # Try to extract username after problem number
    USERNAME=$(echo "$BASENAME" | sed -n "s/.*${PROBLEM_NUM}_\([^_/]*\).*/\1/p" | head -1)
fi

# Final fallback: use directory name or prompt
if [ -z "$USERNAME" ]; then
    if [ -d "$SUBMISSION_INPUT" ]; then
        USERNAME=$(basename "$SUBMISSION_INPUT")
    else
        echo "Error: Could not determine username. Please specify:"
        read -p "Username: " USERNAME
        if [ -z "$USERNAME" ]; then
            echo "Error: Username required"
            exit 1
        fi
    fi
fi

echo "Detected username: $USERNAME"

# Validate required files
if [ ! -f "$TEMP_DIR/policy.json" ]; then
    echo "Error: policy.json not found in submission"
    exit 1
fi

if [ ! -f "$TEMP_DIR/submission.json" ]; then
    echo "Warning: submission.json not found, creating minimal one..."
    python3 << EOF
import json
with open("$TEMP_DIR/submission.json", "w") as f:
    json.dump({
        "problem_id": "$PROBLEM_ID",
        "team_name": "$USERNAME",
        "submission_date": "$(date +%Y-%m-%d)",
        "policy_version": "1.0"
    }, f, indent=2)
EOF
fi

# Create target directory
TARGET_DIR="submissions/$PROBLEM_ID/$USERNAME"
mkdir -p "$TARGET_DIR"

# Copy files
echo "Organizing submission into: $TARGET_DIR"
cp "$TEMP_DIR"/*.json "$TARGET_DIR/" 2>/dev/null || true

# Verify
if [ ! -f "$TARGET_DIR/policy.json" ]; then
    echo "Error: Failed to copy policy.json"
    exit 1
fi

echo ""
echo "=========================================="
echo "Import Complete!"
echo "=========================================="
echo "Submission organized in: $TARGET_DIR"
echo ""
echo "Next step: Process the submission"
echo "  tools/process_submission.sh $TARGET_DIR"
echo ""

