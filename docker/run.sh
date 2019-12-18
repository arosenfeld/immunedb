#!/usr/bin/env bash

echo 'Moving MySQL to Volume'
mkdir -p /share/mysql_data /share/configs
chown mysql:mysql /share/mysql_data
rsync -a --ignore-existing /var/lib/mysql/* /share/mysql_data
service mysql start
echo 'Setting up database'
cat /tmp/setup_users.sql | mysql -u root &> /dev/null
echo 'Starting webserver'
if [ -z "$IMMUNEDB_DAEMON" ]
then
    python3 /apps/immunedb/proxy.py &> /apps/immunedb/proxy.log &
else
    python3 /apps/immunedb/proxy.py
fi
