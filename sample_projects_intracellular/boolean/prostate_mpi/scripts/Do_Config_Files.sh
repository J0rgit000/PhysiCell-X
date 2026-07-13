#!/bin/bash

#if [ "$#" -lt 3 ]; then
#    echo "Uso: $0 <origen_file> <simulations_number> <project_name>"
#    exit 1
#fi

# Variables
ARCHIVO_ORIGEN=$(readlink -f "$1")
MAX_SIMS=${2:-2} #We put 2 simulations as default
SIMULADOR=$3

# Output script in the current directory
FILENAME="$(pwd)/run_drug_simulation.sh"

if [ ! -f "$ARCHIVO_ORIGEN" ]; then
    echo "Error: The file '$ARCHIVO_ORIGEN' doesn't exist."
    exit 1
fi

NOMBRE_BASE="${ARCHIVO_ORIGEN%.*}"
EXTENSION="${ARCHIVO_ORIGEN##*.}"

echo "Preparing '$MAX_SIMS' simulations for '$ARCHIVO_ORIGEN'..."

for i in $(seq 1 "$MAX_SIMS")
do
    NAME_ARCHIVO="${NOMBRE_BASE}_$i"
    ARCHIVO_TEMP="${NOMBRE_BASE}_$i.${EXTENSION}"
    
    cp "$ARCHIVO_ORIGEN" "$ARCHIVO_TEMP"

    # 1. Update internal folder path in the file
    sed -i -E "s#(<folder>output/[^/]*/)([^<]*)#\1\2_$i#g" "$ARCHIVO_TEMP"
    
    # 2. Create output directory logic
    NEW_PATH=$(echo "$NAME_ARCHIVO" | sed 's|^.*/config/|output/|')
    mkdir -p "$NEW_PATH"
    
    # 3. CONVERT TO PARTIAL PATH FOR THE RUN SCRIPT
    # This removes everything up to "config/"
    PARTIAL_PATH="config/${ARCHIVO_TEMP#*/config/}"
    
    # Append the clean path to the runner script
    echo "./$SIMULADOR $PARTIAL_PATH" >> "$FILENAME"
done
