FROM ghcr.io/aiidateam/aiida-core-with-services:latest

USER root

ARG QE_VERSION
ENV QE_VERSION ${QE_VERSION}

# Install aiida-quantumespresso from source code
COPY --from=src . /tmp/aiida-quantumespresso
RUN pip install /tmp/aiida-quantumespresso --no-cache-dir && \
    rm -rf /tmp/aiida-quantumespresso

# Install quantum espresso from conda (the latest version)
RUN mamba install -y -c conda-forge qe=${QE_VERSION} && \
    mamba clean --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${SYSTEM_USER}"

COPY scripts/60-pseudo-code-setup.sh /etc/init/run-before-daemon-start/

# Static example files
RUN mkdir -p /opt/examples
COPY _static/ /opt/examples

USER ${SYSTEM_UID}
