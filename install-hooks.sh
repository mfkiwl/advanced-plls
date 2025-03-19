#!/bin/bash
# This script sets Git to use the custom hooks in the githooks folder.

git config core.hooksPath githooks
echo "Custom hooks installed. Your repository now uses the hooks in the 'githooks' folder."
