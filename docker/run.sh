#!/usr/bin/env bash

echo 'Moving MySQL to Volume'
mkdir -p /share/mysql_data /share/configs
chown mysql:mysql /share/mysql_data
rsync -a --ignore-existing /var/lib/mysql/* /share/mysql_data
service mysql start
echo 'Setting up database'
cat /tmp/setup_users.sql | mysql -u root &> /dev/null
#immunedb_admin create immunedb /root/configs --admin-pass insecure_password
#immunedb_rest /root/configs/immunedb.json &> /var/log/immunedb-rest.log &
