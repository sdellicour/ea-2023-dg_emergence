#!/bin/bash
#SBATCH --job-name=NA_030624_3
#SBATCH --time=119:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=10240
#SBATCH --qos=gpu_prio

module load releases/2022a
module load beagle-lib/4.0.0-GCC-11.3.0-CUDA-11.7.0

java -jar beast_dev_1105_290920.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite NA_030624_3.xml
