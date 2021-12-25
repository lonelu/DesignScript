-----------------------------------
#Follow the manual here: 
Following PyRosetta Tests FAILED: 
http://www.pyrosetta.org/documentation/windows
http://www.pyrosetta.org/dow

#set up pymol:
http://www.pyrosetta.org/pymol_mover-tutorial
#Download pymol with python2.7 supported. https://pymol.org/2/  ->  https://pymol.org/installers/

-----------------------------------

#How to setup sublinux and pyrosetta.

> sudo apt update
> sudo apt install python3-pip
> sudo apt install gcc #(not sure about if this cmd is necessary)

#install miniconda.
> bash Miniconda3-latest-Linux-x86_64.sh

#create and activate conda environmen
> conda env create -f environment_linux.xml
> conda activate env_conda

#install pyrosetta
> tar -vjxf PyRosetta-<version>.tar.bz2
> python3 PyRosetta/setup/setup.py install

> nano test_pyrosetta_install.py
>>from pyrosetta import *
>>init()

> python3 test_pyrosetta_install.py

