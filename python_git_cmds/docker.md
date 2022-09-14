### Dockerfile from https://hub.docker.com/r/230218818/colabfold

FROM nvidia/cuda:11.1.1-devel-ubuntu20.04

COPY ./sources.list /etc/apt/sources.list
COPY ./install_colabbatch_linux.sh /

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys A4B469963BF863CC \
&& apt update \
&& apt -y install gcc curl git wget gnutls-bin \
&& git config --global http.sslVerify false \
&& git config --global http.postBuffer 1048576000 \
&& bash /install_colabbatch_linux.sh

ENV PATH="/colabfold_batch/conda/bin:/colabfold_batch/conda/condabin:/colabfold_batch/bin:$PATH"
RUN echo "source activate /colabfold_batch/colabfold-conda" > ~/.bashrc

CMD /bin/bash

### Build by pull in gpu
sudo docker pull 230218818/colabfold

## How to run it in gpu.
sudo docker run -d -v /home/gpu/Lei/DesignData/ntf2/:/home/colabfold/ --rm --gpus all 230218818/colabfold /bin/bash -c "source activate /colabfold_batch/colabfold-conda && colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax /home/colabfold/input/all_seq.fasta /home/colabfold/output2/"

sudo docker run -it -v /home/gpu/Lei/DesignData/ntf2/:/home/colabfold/ 230218818/colabfold /bin/bash 
colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax /home/colabfold/input/all_seq.fasta /home/colabfold/output2/ --cpu


## Build from Docker file in gpu
FROM nvidia/cuda:11.1.1-devel-ubuntu20.04

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys A4B469963BF863CC \
&& apt update \
&& apt -y install gcc curl git wget gnutls-bin \
&& git config --global http.sslVerify false \
&& git config --global http.postBuffer 1048576000 \
&& wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
&& bash /install_colabbatch_linux.sh

ENV PATH="/colabfold_batch/conda/bin:/colabfold_batch/conda/condabin:/colabfold_batch/bin:$PATH"
RUN echo "source activate /colabfold_batch/colabfold-conda" > ~/.bashrc

CMD /bin/bash

## Run
sudo docker run -d -v /home/gpu/Lei/DesignData/ntf2/:/home/ --rm --gpus all 9c361110565c /bin/bash -c "source activate /colabfold_batch/colabfold-conda && colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax /home/colabfold/input/all_seq.fasta /home/colabfold/output3/"

## Run interactive. The following works but get Value Error."ValueError: Minimization failed after 100 attempts."
sudo docker run -it -v /home/gpu/Lei/DesignData/ntf2/:/home/ 8c37c2bb8a18 /bin/bash
colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax /home/input/all_seq.fasta /home/colabfold/output2/ --cpu


### On wynton build by pull
singularity build docker_colabfold.sif docker://230218818/colabfold

## Run interactive
singularity run rocker_r-base.sif

colabfold_batch /wynton/home/degradolab/lonelu/DesignData/NTF2/logx_tts_r3_colabfold/all_seq.fasta /wynton/home/degradolab/lonelu/DesignData/NTF2/logx_tts_r3_colabfold/outdir/ --cpu

## Got the following error. Can we download it when we build it.
FileNotFoundError: [Errno 2] No such file or directory: '//colabfold_batch/colabfold/params'

#
2022-06-09 14:09:37,676 Running colabfold 1.3.0 (03b9a6ad72c1ea90173f55d252ff25dd4420906d)
Traceback (most recent call last):
  File "/colabfold_batch/colabfold-conda/lib/python3.7/pathlib.py", line 1273, in mkdir
    self._accessor.mkdir(self, mode)
FileNotFoundError: [Errno 2] No such file or directory: '//colabfold_batch/colabfold/params             '

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "//colabfold_batch/colabfold-conda/bin/colabfold_batch", line 8, in <module>
    sys.exit(main())
  File "/colabfold_batch/colabfold-conda/lib/python3.7/site-packages/colabfold/batch.py", l             ine 1679, in main
    download_alphafold_params(model_type, data_dir)
  File "/colabfold_batch/colabfold-conda/lib/python3.7/site-packages/colabfold/download.py"             , line 36, in download_alphafold_params
    params_dir.mkdir(parents=True, exist_ok=True)
  File "/colabfold_batch/colabfold-conda/lib/python3.7/pathlib.py", line 1277, in mkdir
    self.parent.mkdir(parents=True, exist_ok=True)
  File "/colabfold_batch/colabfold-conda/lib/python3.7/pathlib.py", line 1273, in mkdir
    self._accessor.mkdir(self, mode)
OSError: [Errno 30] Read-only file system: '//colabfold_batch/colabfold'

#


Downloading alphafold2 weights to //colabfold_batch/colabfold:   0%| | 10.0k/3.47G [00:00<1
Traceback (most recent call last):
  File "//colabfold_batch/colabfold-conda/bin/colabfold_batch", line 8, in <module>
    sys.exit(main())
  File "/colabfold_batch/colabfold-conda/lib/python3.7/site-packages/colabfold/batch.py", l             ine 1681, in main
    download_alphafold_params(model_type, data_dir)
  File "/colabfold_batch/colabfold-conda/lib/python3.7/site-packages/colabfold/download.py"             , line 47, in download_alphafold_params
    file.extractall(path=params_dir)
  File "/colabfold_batch/colabfold-conda/lib/python3.7/tarfile.py", line 2002, in extractal             l
    numeric_owner=numeric_owner)
  File "/colabfold_batch/colabfold-conda/lib/python3.7/tarfile.py", line 2044, in extract
    numeric_owner=numeric_owner)
  File "/colabfold_batch/colabfold-conda/lib/python3.7/tarfile.py", line 2114, in _extract_             member
    self.makefile(tarinfo, targetpath)
  File "/colabfold_batch/colabfold-conda/lib/python3.7/tarfile.py", line 2155, in makefile
    with bltn_open(targetpath, "wb") as target:
OSError: [Errno 30] Read-only file system: '//colabfold_batch/colabfold/params/./params_mod             el_2.npz'

