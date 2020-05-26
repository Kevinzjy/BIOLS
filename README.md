# BIOLS
Useful tools on BIOLS HPC

# Python programming tools

```python
import sys
sys.path.append('/home/zhangjy/git/BIOLS')
import py3biotools
```

# Bioinformatic Tools
- bin/FastqToFasta : Convert FastQ to FastA format
- bin/fastcount    : count total base and N50 of input fasta/fastq

# Userful linux Tools:
- bin/FindHulk     : Find large space-consuming files in directory

# PBS server tools

## scripts
- bin/job-db-backup.cron : cron job for backing up PBS server accounting log
- bin/pbs_accounting     : auto accouting for each PBS user, python2 only
- bin/pestat             : node status monitor for PBS
- bin/pgrep_watch        : send notification when a process finished
- bin/pushover           : push notification to my PUSHOVER account
- bin/qdel_all           : delete all jobs of a user
- bin/rpm_unpack         : unpack rpm packages
- bin/server_status      : check status of PBS server when login (add to .bash_profile)
- bin/user_space         : account space usage for all user (home directory only), python2 only

## configs
- HPC/ban_ip.sh          : ban certain ip using firewall-cmd/fail2ban
- HPC/new_user.sh        : batch add newuser to add nodes using pdsh
- HPC/torque.cfg         : Torque config

Remember to add PDSH environment variable for using ssh

```bash
export PDSH_RCMD_TYPE=ssh
```

