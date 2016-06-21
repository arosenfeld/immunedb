#!/usr/bin/env bash
set -e
function setup() {
    immunedb_admin create test_db . --admin-pass "$DB_ADMIN_PASS"
}

function teardown() {
    immunedb_admin delete test_db.json --delete-user --admin-pass "$DB_ADMIN_PASS"
    rm test_db.json
}

if [ -z "$LL_PATH" ]
then
    echo 'LL_PATH must be set'
    exit
fi

if [ -z "$NO_TEARDOWN" ]
then
    trap teardown 0
else
    echo 'Not tearing down since NO_TEARDOWN is set'
fi

setup
coverage erase
coverage run --source=immunedb -p -m nose -s tests/tests_parser.py
coverage run --source=immunedb -p -m nose -s tests/tests_pipeline.py
coverage run --source=immunedb -p -m nose -s tests/tests_import.py
