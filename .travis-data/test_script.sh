#!/bin/bash

# Be verbose, and stop with error as soon there's one
set -ev

case "$TEST_TYPE" in
    docs)
        # Compile the docs (HTML format); -W to convert warnings in errors, 
        # -n to warn about all missing references
        SPHINXOPTS="-nW" make -C docs html
        ;;
    tests)
        verdi -p test_$TEST_AIIDA_BACKEND devel tests db.quantumespresso

        # Run the daemon tests using docker
        verdi -p $TEST_AIIDA_BACKEND run ${TRAVIS_BUILD_DIR}/.travis-data/test_pw_with_daemon.py
        ;;
    pre-commit)
        pre-commit run --all-files || ( git status --short ; git diff ; exit 1 )
        ;;
esac
