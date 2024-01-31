# Environment
All scripts were run under Ubuntu 22.04.2 LTS Linux

## Dependencies
- Python 3.11
- CUDA Toolkit 11.5 (https://developer.nvidia.com/cuda-11-5-0-download-archive)
- Boost 1.84.0 (https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.bz2)
- Vina-GPU 2.0 (https://github.com/DeltaGroupNJUPT/QuickVina2-GPU)

# Setup
Install python dependencies by running
```bash
pip install --trusted-host pypi.python.org -r src/requirements.txt
```
# Run high-throughput virtual screening
Run screening by running
```bash
./run_screening.sh <PATH_TO_BOOST_LIBRARY> <PATH_TO_CUDA_LIBRARIES> <PATH_TO_VINA_GPU>
```
