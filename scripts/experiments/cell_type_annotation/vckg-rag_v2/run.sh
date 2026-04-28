#!/bin/bash

# Example GATHER run. Override paths and parameters from the command line
# by calling run_inference.py directly.

python run_inference.py \
    --top-k-genes 50 \
    --top-k-celltypes 15 \
    --max-hops 2 \
    --max-concurrent 100 
