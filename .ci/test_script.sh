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
        pytest -v ${TRAVIS_BUILD_DIR}/tests
        ;;
    pre-commit)
        pre-commit run --all-files || ( git status --short ; git diff ; exit 1 )
        ;;
esac
