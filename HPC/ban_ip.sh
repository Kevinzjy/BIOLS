IP=$1
firewall-cmd --permanent --add-rich-rule="rule family='ipv4' source address='$IP' reject"

fail2ban-client status sshd
fail2ban-client set sshd unbanip 172.16.100.1
