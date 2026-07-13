#!/bin/bash

# Get script location to find Prepare_Jobs.sh
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORKER_SCRIPT="$SCRIPT_DIR/Do_Config_Files.sh"

#if [ "$#" -lt 3 ]; then
#    echo "Uso: $0 <directory> <simulations_number> <project_name>"
#    echo "Not enough parameters for Run_Basic_Config_Files.sh"
#    exit 1
#fi

# Convert the input directory to a full absolute path
TARGET_DIR=$(readlink -f "$1")
SIMNUMBER=$2
PROJECT=$3

if [ ! -d "$TARGET_DIR" ]; then
  echo "Error: Directory '$TARGET_DIR' does not exist."
  exit 1
fi

for file in "$TARGET_DIR"/*; do
  if [ -f "$file" ]; then
    #echo "Sending $file to worker..."
    bash "$WORKER_SCRIPT" "$file" "$SIMNUMBER" "$PROJECT"
  fi
done
