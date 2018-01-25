#!/bin/sh

echo "Running GCAM"

OUTPUT="/results"
INPUT="/input.txt"
DATABASE="/resource"

if [ -z "$1" ]; then
    echo "No additional arguments for run...."
    python gcamrun genebased "-o" "$OUTPUT" "-p" "$INPUT" "-d" "$DATABASE"
elif [ -z "$2" ]; then
    echo "There is one arguments: $1"
    python gcamrun genebased "-o" "$OUTPUT" "-p" "$INPUT" "-d" "$DATABASE" "$1"
else
    echo "There are two arguments: $1 $2"
    python gcamrun genebased "-o" "$OUTPUT" "-p" "$INPUT" "-d" "$DATABASE" "$1" "$2"
fi
