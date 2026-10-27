To run STRONG follow these steps:

0. Get an account on UoWs high performance cluster Avon.

1. You will need to move files onto and from the cluster. You can do this by running "sftp <username>@myfiles.warwick.ac.uk"

2. Follow the instructions here:
https://github.com/OSS-Lab/STRONG_pipeline_avon

3. To execute command files, such as "config yaml", you will need to use the sbatch files found in this folder. 

4. Running a strong workflow will require you to transfer files over to the server and also have them in a format where you have only two sets of files for each sample (assuming paired end sequencing). You place the files in a folder structure which needs to be described in the "config.yaml" file - see example. Also, if you have more than 2 files per sample, you can concatenate them into 2 files by using "mergeReads.sh" and "runMergeReeads.sbatch"

5. Jobs on our server have a 48h time limit. It is possible/likely that a Strong workflow will not finish in that time. In that case, you can simply restart the job and things will continue from where it was left. However, you should first run "runStrongUnlock.sbatch" to unlock any output folders.