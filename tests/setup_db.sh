#!/usr/bin/env bash
echo "innodb_large_prefix = ON" >> /etc/mysql/my.cnf
echo "innodb_file_format = Barracuda" >> /etc/mysql/my.cnf
sldb_admin create root test_db . --admin-pass ''
