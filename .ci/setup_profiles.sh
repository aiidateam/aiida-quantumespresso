#!/usr/bin/env bash
set -ev

if [[ "$TEST_TYPE" != "pre-commit" ]]
    # no setup at all required for pre-commit to run
then
    reentry scan

    # Create the main and test database
    psql -h localhost -c "CREATE DATABASE $TEST_AIIDA_BACKEND;" -U postgres -w
    psql -h localhost -c "CREATE DATABASE test_$TEST_AIIDA_BACKEND;" -U postgres -w

    # Setup the main profile
    verdi setup  --profile $TEST_AIIDA_BACKEND \
        --email='aiida@localhost' --first-name='AiiDA' --last-name='test' --institution='AiiDA Team' \
        --db-engine='postgresql_psycopg2' --db-backend="${TEST_AIIDA_BACKEND}" --db-host='localhost' --db-port=5432 \
        --db-name="$TEST_AIIDA_BACKEND" --db-username=postgres --db-password='' \
        --repository="/tmp/repository_${TEST_AIIDA_BACKEND}/" --non-interactive

    # Setup the test profile
    verdi setup --profile test_$TEST_AIIDA_BACKEND \
        --email='aiida@localhost' --first-name='AiiDA' --last-name='test' --institution='AiiDA Team' \
        --db-engine='postgresql_psycopg2' --db-backend="${TEST_AIIDA_BACKEND}" --db-host='localhost' --db-port=5432 \
        --db-name="test_$TEST_AIIDA_BACKEND" --db-username=postgres --db-password='' \
        --repository="/tmp/test_repository_test_${TEST_AIIDA_BACKEND}/" --non-interactive

    verdi profile setdefault $TEST_AIIDA_BACKEND
    verdi config runner.poll.interval 0
fi
