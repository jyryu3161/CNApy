#!/bin/bash
#
# CNApy launcher script
# Usage: ./run.sh [project_path] [scenario_path]
#

# Move to the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Check the virtual environment and run via uv
echo "🧬 Starting CNApy..."
echo "   Project directory: $SCRIPT_DIR"

# Disable Qt debug logging
export QT_LOGGING_RULES="*.debug=false"

# Check whether uv is installed
if command -v uv &> /dev/null; then
    echo "   Running via uv..."
    uv run python cnapy.py "$@"
else
    # Fall back to .venv if uv is not available
    if [ -d ".venv" ]; then
        echo "   Running via virtual environment (.venv)..."
        source .venv/bin/activate
        python cnapy.py "$@"
    else
        echo "❌ Error: uv is not installed and no .venv virtual environment was found."
        echo "   Please run 'uv sync' first to set up the environment."
        exit 1
    fi
fi
