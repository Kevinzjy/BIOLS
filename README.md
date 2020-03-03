# BIOLS
Useful tools on BIOLS HPC

# Python programming tools
- python/utils.py : useful tools
- python/bioinfo.py : python bioinformatic tools

# Bioinformatic Tools
- FastqToFasta : Convert FastQ to FastA format
- FindHulk     : Find large space-consuming files in directory
- fastcount    : count total base and N50 of input fasta/fastq

# PBS server maintenance
- job-db-backup.cron : cron job for backing up PBS server accounting log
- pbs_accounting     : auto accouting for each PBS user, python2 only
- pestat             : node status monitor for PBS
- pgrep_watch        : send notification when a process finished
- pushover           : push notification to my PUSHOVER account
- qdel_all           : delete all jobs of a user
- rpm_unpack         : unpack rpm packages
- server_status      : check status of PBS server when login (add to .bash_profile)
- user_space         : account space usage for all user (home directory only), python2 only
