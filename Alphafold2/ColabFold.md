### To install colabfold in local please follow:
https://github.com/YoshitakaMo/localcolabfold
or 
https://github.com/sokrypton/ColabFold (Prefer this way of installation. However, the setup_databases.sh is not working.)

# Install in wynton
>>> ssh dev1 or ssh gpudev1

Note that wynton default gcc is 4.8.5, please use devtoolset: https://wynton.ucsf.edu/hpc/software/scl.html#developer-toolset-scls
>>> scl enable devtoolset-7 bash 
>>> gcc --version
Check the gcc version should be different now.

>>> wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
>>> bash install_colabbatch_linux.sh

>>> export PATH="/wynton/home/degradolab/lonelu/GitHub_Design/colabfold_batch/bin:$PATH"
Paste this in .bashrc

# Install in wynton


# Run prediction.
environment location: /wynton/home/degradolab/lonelu/GitHub_Design/colabfold_batch/colabfold-conda

>>> colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax inputfile outputdir/ 


# Note: 
It seems like the cuda version is still not correct.
'''
    You may not need to update to CUDA 11.1; cherry-picking the ptxas binary is often sufficient.
    2022-08-05 16:37:14.832263: W external/org_tensorflow/tensorflow/stream_executor/gpu/asm_compiler.cc:111] *** WARNING *** You are using ptxas 9.1.121, which is older than 11.1. ptxas before 11.1 is known to miscompile XLA code, leading to incorrect results or invalid-address errors
'''

Overall, not fully successful.
