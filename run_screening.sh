#!/bin/bash

usage() {
    echo "Usage: ./run_screening.sh <PATH_TO_BOOST> <PATH_TO_CUDA_LIBRARIES> <PATH_TO_VINA_GPU>"
    exit 1
}
if [ "$#" -ne 3 ]; then
    usage
fi

OLD_PATH=${PATH}
OLD_LD=${PATH}

export PATH=${PATH}:${$1}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${$2}

python3 src/main.py $3

export PATH=${OLD_PATH}
export LD_LIBRARY_PATH=${OLD_LD}
