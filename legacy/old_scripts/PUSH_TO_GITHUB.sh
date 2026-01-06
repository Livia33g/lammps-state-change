#!/bin/bash
# Script to push this repository to GitHub
# Usage: ./PUSH_TO_GITHUB.sh YOUR_USERNAME REPO_NAME

if [ $# -ne 2 ]; then
    echo "Usage: $0 YOUR_USERNAME REPO_NAME"
    echo "Example: $0 Livia33g lammps-state-change"
    exit 1
fi

USERNAME=$1
REPO_NAME=$2

echo "Adding GitHub remote..."
git remote add origin https://github.com/${USERNAME}/${REPO_NAME}.git

echo "Setting branch to main..."
git branch -M main

echo "Pushing to GitHub..."
git push -u origin main

echo "✅ Repository pushed to GitHub!"
echo "View it at: https://github.com/${USERNAME}/${REPO_NAME}"

