#!/usr/bin/env bash
function setup() {
    sldb_admin create root test_db . --admin-pass "$DB_ADMIN_PASS"
}

function teardown() {
    sldb_admin delete root test_db.json --delete-user --admin-pass "$DB_ADMIN_PASS"
    rm test_db.json
}

setup
nosetests -s tests/tests_pipeline.py
#nosetests -s tests/tests_import.py

if [ -z "$NO_TEARDOWN" ]
then
    teardown
else
    echo 'Not tearing down since NO_TEARDOWN is set'
fi
