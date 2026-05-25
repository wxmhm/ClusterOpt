#!/bin/bash
# ============================================================================
# ClusterOpt Project Backup Script
# Dynamically discovers files via git - no hardcoded file lists.
# New .cpp/.h/data files are automatically included.
# ============================================================================
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
BACKUP_NAME="ClusterOpt_$(date +%Y%m%d_%H%M%S).zip"

# ------------------------------------------------------------------
# Exclusion patterns (edit these to adjust what gets excluded)
# These are grep -v patterns applied to the file list.
# ------------------------------------------------------------------
EXCLUDE_PATTERNS=(
    '^Cluster/results/'
    '^Cluster/Debug/'
    '^Cluster/Release/'
    '^Cluster/ipch/'
    '^Debug/'
    '^Release/'
    '^x64/'
    '^ipch/'
    '^\.vs/'
    '\.sdf$'
    '\.opendb$'
    '\.VC\.db$'
)

cd "$PROJECT_ROOT"

# Build grep filter from patterns
GREP_FILTER=""
for pat in "${EXCLUDE_PATTERNS[@]}"; do
    if [ -z "$GREP_FILTER" ]; then
        GREP_FILTER="grep -v '$pat'"
    else
        GREP_FILTER="$GREP_FILTER | grep -v '$pat'"
    fi
done

# Collect all files known to git (tracked + untracked-but-not-ignored)
echo "Discovering files..."
FILE_LIST=$( {
    git ls-files --cached
    git ls-files --others --exclude-standard
} | sort -u | eval "$GREP_FILTER")

FILE_COUNT=$(echo "$FILE_LIST" | grep -c . || true)
echo "Found $FILE_COUNT files to backup."

# Create archive
echo "Creating $BACKUP_NAME ..."
echo "$FILE_LIST" | zip -q -@ "$BACKUP_NAME"

SIZE=$(du -h "$BACKUP_NAME" | cut -f1)
echo "Done: $BACKUP_NAME ($SIZE)"
echo ""
echo "Files included in backup:"
echo "$FILE_LIST" | head -20
if [ "$FILE_COUNT" -gt 20 ]; then
    echo "... and $((FILE_COUNT - 20)) more files"
fi
echo ""
echo "Restore with: unzip $BACKUP_NAME -d <target_dir>"
