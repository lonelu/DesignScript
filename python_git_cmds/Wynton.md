# If there is a problem log in lonelu@log2.wynton.ucsf.edu
# Open the config file: C:\Users\lulei_000\.ssh.\known_hosts Remove the related key.

# Submit a job. Example:
> cd ~/ab_initio_example
# change the filepath in abinitio_G8.job
> qsub abinitio_G8.job
# check the output

/wynton/home/degradolab/smann2/programs/rosetta_bin_linux_2018.33.60351_bundle/main/source/bin/combine_silent.static.linuxgccrelease -in:file:silent orig_abinitio_silent_*.out -in:file:silent_struct_type binary -out:file:silent orig_abinitio_silent_combined.out -out:file:silent_struct_type binary



