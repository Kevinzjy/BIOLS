new_uid=2000
new_user="user_id"
new_home="user_home"
new_group="user_group"

pdsh -g compute_nodes -l root "useradd -g ${new_group} -u ${new_uid} -d ${new_home} ${new_user}"

