#!/bin/bash

# Be verbose, and stop with error as soon there's one
set -ev

if [[ "$TEST_TYPE" == "tests" ]]
then
    THEKEY=`cat ${HOME}/.ssh/id_rsa.pub`
    echo 'AUTHORIZED_KEY='"$THEKEY" > ${TRAVIS_BUILD_DIR}/torquessh.env
    docker build -t torquessh ${TRAVIS_BUILD_DIR}/.ci/torquessh-qe
    # Run it in the background, mapping port 22 of the container
    # to port 10022 outside, and passing the environment variable
    docker run -d --privileged -p=10022:22 --env-file ${TRAVIS_BUILD_DIR}/torquessh.env torquessh
    # Docker ps to see what is going on
    echo "Running docker ps to see if the 'torquessh' docker image is up..."
    docker ps
    # Wait for SSH to be up
    ${TRAVIS_BUILD_DIR}/.ci/wait-for-it.sh localhost:10022 -t 0

fi
