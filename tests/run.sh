#!/usr/bin/env bash
function setup() {
    sldb_admin create root test_db . --admin-pass "$DB_ADMIN_PASS"
}

function teardown() {
    sldb_admin delete root test_db.json --delete-user --admin-pass "$DB_ADMIN_PASS"
    rm test_db.json
}

setup
nosetests -s
teardown
