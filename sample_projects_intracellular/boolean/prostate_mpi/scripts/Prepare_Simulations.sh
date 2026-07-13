#!/bin/bash

# Get the absolute path
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORKER_SCRIPT="$SCRIPT_DIR/Run_Basic_Config_Files.sh"

#I need this two parameters to pass them to the final script.
if [ "$#" -lt 2 ]; then
    echo "Parameters to give: $0 <simulations_number> <project_name>"
    echo "Ex: $0 10 prostate"
    exit 1
fi

SIM_NUM=$1
PROJECT_NAME=$2

echo "Starting sequential processing from: $SCRIPT_DIR"

# Process folders (using relative paths from where you run the script)
echo "[$(date +%T)] --- Processing Single folder ---" #Single path hard coded as this is think for MaBoSS.
bash "$WORKER_SCRIPT" "config/single" "$SIM_NUM" "$PROJECT_NAME"

echo "[$(date +%T)] --- Processing Double folder ---" 
bash "$WORKER_SCRIPT" "config/double" "$SIM_NUM" "$PROJECT_NAME"

echo "Finished all tasks, all simulations set up clap clap."
