#!/usr/bin/env bash
sldb_admin create root test_db . --admin-pass ''

if [ `mysql -u root -e 'show databases' | grep test_db | wc -l` != 1 ]
then
    exit -1
fi

