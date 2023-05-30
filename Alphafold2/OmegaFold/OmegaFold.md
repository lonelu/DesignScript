So OmegaFold works on Wynton, where there are GPUs with enough memory to store the parameters.  To get it to work there, I recommend taking the following steps:
1. git clone the OmegaFold repository from GitHub (https://github.com/HeliXonProtein/OmegaFold) (must be on a dev node to use git)
2. Create a new conda environment for OmegaFold: conda create -n env_omegafold python=3.10.4 (I think you'll have some freedom to choose your Python version here)
3. Install the requirements using pip install -r requirements.txt (once you're in the OmegaFold directory you cloned from GitHub)
4. ssh dt1 and scp -r gpu@10.36.173.60:/home/gpu/.cache/omegafold_ckpt ~/.cache to get the model parameters, then exit  to get back to the dev node
5. Add the following line to OmegaFold/main.py at Line 30: torch.cuda.set_per_process_memory_fraction(0.8)
6. Download the submit script attached to this Slack message to your OmegaFold directory on Wynton and run OmegaFold with the following command: qsub -q gpu.q submit_main.sh <PATH TO YOUR FASTA> <PATH TO YOUR OUTPUT DIRECTORY>.  This will run just as main.py would run locally if the dev node had a GPU.

After step3 check the github readme. install. Always note that the program can be updated.