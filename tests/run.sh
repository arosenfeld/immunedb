#!/usr/bin/env bash
set -e
function setup() {
    immunedb_admin create test_db . --admin-pass "$DB_ADMIN_PASS"
}

function teardown() {
    trap '' INT TERM
    immunedb_admin delete test_db.json --delete-user --admin-pass "$DB_ADMIN_PASS"
    rm test_db.json
}

if [ -z "$NO_TEARDOWN" ]
then
    trap teardown 0
else
    echo 'Not tearing down since NO_TEARDOWN is set'
fi

setup
coverage erase
coverage run --source=immunedb -p -m nose tests/tests_parser.py
coverage run --source=immunedb -p -m nose tests/tests_import.py
coverage run --source=immunedb -p -m nose tests/tests_pipeline.py
coverage run --source=immunedb -p -m nose tests/run_server.py &
PID=$!
sleep 5
coverage run --source=immunedb -p -m nose tests/tests_api.py
coverage run --source=immunedb -p -m nose tests/tests_export.py
coverage run --source=immunedb -p -m nose tests/tests_clone_import.py
