#!/bin/bash
# Script to retrieve relevant distance data from .xtc trajectories
# Author: Boris N. Schüpp
# Check if both arguments are provided
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
  echo "Usage: $0 <Trajectory directory> <Structure directory> <End to end distance ouput directory>"
  exit 1
fi

TRAJECTORY_DIRECTORY=$1
STRUCTURE_DIRECTORY=$2
OUTPUT_DIRECTORY=$3

# Check if the trajectory directory is a valid directory
if [ ! -d "$TRAJECTORY_DIRECTORY" ]; then
  echo "Error: '$TRAJECTORY_DIRECTORY' is not a valid directory."
  exit 1
fi

# Check if the structure directory is a valid directory
if [ ! -d "$STRUCTURE_DIRECTORY" ]; then
  echo "Error: '$STRUCTURE_DIRECTORY' is not a valid directory."
  exit 1
fi

# Iterate over files in the trajectory directory
echo "Processing files in the Trajectory directory: $TRAJECTORY_DIRECTORY"
for file in "$TRAJECTORY_DIRECTORY"/*.xtc; do
	echo "Processing $file"
	filename=$(basename "$file" .xtc)
	
	[[ ! "$filename" =~ ^[01] ]] && continue

	if [[ "$file" == 0* ]]; then
		python DetermineIndicesForceVariation.py "0_GTT100.gro" "${STRUCTURE_DIRECTORY}"
	else
		python DetermineIndicesForceVariation.py "1_GTT100.gro" "${STRUCTURE_DIRECTORY}"
	fi

	mapfile -t arr < <(seq 0 $(( $(wc -l <TempBaseDist.ndx) / 2 - 1 )))
	mapfile -t arr2 < <(seq 0 $(( $(wc -l <TempBackboneDist.ndx) / 2 - 1 )))
	# Extract the distances using the created indices
	gmx distance -f "$file" -s "${STRUCTURE_DIRECTORY}"/"${filename}".tpr -n TempEtE.ndx -oall TempEndToEndDistances.xvg -select 0 -nopbc
	gmx distance -f "$file" -s "${STRUCTURE_DIRECTORY}"/"${filename}".tpr -n TempBaseDist.ndx -oall TempBaseDistances.xvg -select "${arr[@]}"
	gmx distance -f "$file" -s "${STRUCTURE_DIRECTORY}"/"${filename}".tpr -n TempBackboneDist.ndx -oall TempBackboneDistances.xvg -select "${arr2[@]}" -b 20000 > TempBackboneDistances.txt
	# Process the obtained distances
	python ProcessData.py "$OUTPUT_DIRECTORY" "$filename" "105"
	# Cleanup temporary files
	rm Temp*
done 

# Combine results to single files containing info on all sequence, nick, run pairs. Cleanup the single files.
python CombineResults.py "$OUTPUT_DIRECTORY"
find "$OUTPUT_DIRECTORY" -type f ! -name '*all*' -exec rm {} \;
