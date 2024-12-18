#!/bin/bash

SELECTED=$1

if [ "$SELECTED" = true ]; then
    echo "Selected transcripts"
else
    echo "All transcripts"
fi