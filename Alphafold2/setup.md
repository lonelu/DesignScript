# Try to setup alphafold2 in wynton.
https://www.novopro.cn/articles/202107291210.html

https://github.com/kalininalab/alphafold_non_docker
https://github.com/kuixu/alphafold

# I followed the kuixu's solution. After running the ./install_on_local.sh. I used pip/conda to install several packages.

Finally, I got this.

    (af2) [lonelu@dev1 alphafold]$ python3 run_alphafold_test.py 
    Running tests under Python 3.8.12: /wynton/home/degradolab/lonelu/miniconda3/envs/af2/bin/python3
    2022-01-03 22:07:23.652308: W tensorflow/stream_executor/platform/default/dso_loader.cc:64] Could not load dynamic library 'libcuda.so.1'; dlerror: libcuda.so.1: cannot open shared object file: No such file or directory; LD_LIBRARY_PATH: /wynton/home/degradolab/lonelu/miniconda3/envs/af2/lib
    2022-01-03 22:07:23.652432: W tensorflow/stream_executor/cuda/cuda_driver.cc:269] failed call to cuInit: UNKNOWN ERROR (303)
    [ RUN      ] RunAlphafoldTest.test_end_to_end
    I0103 22:07:23.667741 139853966219072 run_alphafold.py:185] Running model model1
    I0103 22:07:23.668509 139853966219072 run_alphafold.py:195] Total JAX model model1 predict time (includes compilation time, see --benchmark): 0?
    I0103 22:07:23.701562 139853966219072 run_alphafold.py:244] Final timings for test: {'features': 0.00016450881958007812, 'process_features_model1': 8.559226989746094e-05, 'predict_and_compile_model1': 5.173683166503906e-05, 'relax_model1': 0.00010466575622558594}
    [       OK ] RunAlphafoldTest.test_end_to_end
    ----------------------------------------------------------------------
    Ran 1 test in 0.038s

    OK



