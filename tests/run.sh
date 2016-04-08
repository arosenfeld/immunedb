#!/usr/bin/env bash
set -e
function setup() {
    airrdb_admin create root test_db . --admin-pass "$DB_ADMIN_PASS"
}

function teardown() {
    airrdb_admin delete root test_db.json --delete-user --admin-pass "$DB_ADMIN_PASS"
    rm test_db.json
}

if [ -z "$NO_TEARDOWN" ]
then
    trap teardown 0
else
    echo 'Not tearing down since NO_TEARDOWN is set'
fi

setup
nosetests -s tests/tests_pipeline.py
nosetests -s tests/tests_import.py
