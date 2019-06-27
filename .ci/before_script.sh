#!/bin/bash

# Be verbose, and stop with error as soon there's one
set -ev

if [[ "$TEST_TYPE" == "tests" ]]
then
    # Copy test data from the repository to the installed package of the virtual environment
    cp -r ${TRAVIS_BUILD_DIR}/tests ${VIRTUAL_ENV}/lib/python${TRAVIS_PYTHON_VERSION}/site-packages/

    # Refresh the entrypoint cache through reentry
    reentry scan

    # Start the daemon for the correct profile
    verdi -p $TEST_AIIDA_BACKEND daemon start
    verdi -p $TEST_AIIDA_BACKEND daemon status

    # Setup the torquessh computer - this one is custom, using torque scheduler
    verdi -p $TEST_AIIDA_BACKEND computer setup --non-interactive --label=torquessh --hostname=localhost --transport=ssh --scheduler=torque --mpiprocs-per-machine=1 --prepend-text="" --append-text="" --work-dir='/scratch/{username}/aiida_run' --mpirun-command='mpirun -np {tot_num_mpiprocs}'

    # Configure the torquessh computer - this one is custom, using port 22
    verdi -p $TEST_AIIDA_BACKEND computer configure ssh torquessh --non-interactive --safe-interval 0 --username=app --port=10022 --key-filename=~/.ssh/id_rsa --timeout=60 --compress --gss-host=localhost --load-system-host-keys --key-policy=RejectPolicy

    # Configure the 'add' code inside torquessh, which is only required for the integrations test on Jenkins
    verdi -p $TEST_AIIDA_BACKEND code setup -n -L qe-pw \
        -D "Quantum ESPRESSO pw.x code" --on-computer -P quantumespresso.pw \
        -Y torquessh --remote-abs-path=/usr/bin/pw.x

    # Make sure that the torquessh (localhost:10022) key is hashed
    # in the known_hosts file
    echo "'ssh-keyscan -p 10022 -t rsa localhost' output:"
    ssh-keyscan -p 10022 -t rsa localhost > /tmp/localhost10022key.txt
    cat /tmp/localhost10022key.txt

    # Patch for OpenSSH 6, that does not write the port number in the
    # known_hosts file. OpenSSH 7 would work, instead
    if grep -e '^localhost' /tmp/localhost10022key.txt > /dev/null 2>&1 ; then cat /tmp/localhost10022key.txt | sed 's/^localhost/[localhost]:10022/' >> ${HOME}/.ssh/known_hosts ; else  cat /tmp/localhost10022key.txt >> ${HOME}/.ssh/known_hosts; fi

    echo "Content of the known_hosts file:"
    cat ${HOME}/.ssh/known_hosts

fi
