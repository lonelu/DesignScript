# One day my Windows Subsystem started failing with:

The Windows Subsystem for Linux instance has terminated.
The solution I found here was simple enough:

wsl.exe --shutdown
wsl.exe


# Lose access to google drive

sudo mount -t drvfs G: /mnt/g