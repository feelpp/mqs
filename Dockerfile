FROM feelpp/feelpp-toolboxes:feature-biotsavart-v0.108.0-beta.1-ubuntu-20.04

USER root

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# This Dockerfile adds a non-root user with sudo access. Use the "remoteUser"
# property in devcontainer.json to use it. On Linux, the container user's GID/UIDs
# will be updated to match your local UID/GID (when using the dockerFile property).
# See https://aka.ms/vscode-remote/containers/non-root-user for details.
ARG USERNAME=vscode
ARG USER_UID=1001
ARG USER_GID=$USER_UID

ARG BRANCH=biotsavart-bcs

# install Feelpp from BinTray Debian repository
RUN apt update && \
    apt install -y lsb-release && \
    echo "deb http://archive.ubuntu.com/ubuntu $(lsb_release -cs)-proposed main restricted"  |  tee -a /etc/apt/sources.list && \
    echo "deb http://archive.ubuntu.com/ubuntu $(lsb_release -cs)-proposed universe"  |  tee -a /etc/apt/sources.list && \
    echo "deb http://archive.ubuntu.com/ubuntu $(lsb_release -cs)-proposed multiverse"  |  tee -a /etc/apt/sources.list && \
    apt update

# Configure apt and install packages
RUN apt-get update \
    && apt-get -y install --no-install-recommends apt-utils dialog 2>&1 \
    && apt-get -y install build-essential cmake cppcheck valgrind libcurl4-openssl-dev libgsl-dev python3 python3-dev python3-setuptools python3-sympy \
	       libboost1.67-all-dev \
           libcln-dev\
           petsc-dev\
           slepc-dev\
           libhdf5-openmpi-dev\
           libnlopt-dev\
           libgsl-dev\
           libnetcdf-dev libgl2ps-dev libglu1-mesa-dev libsm-dev libxt-dev\
           libfftw3-mpi-dev\
	       libxml2-dev\
	       libgmsh-dev\
	       libtbb-dev\
	       libann-dev libglpk-dev\
           libbz2-dev\
           libbson-dev\
           libmongoc-dev\
           libmongoclient-dev           \
           libglew-dev
    
# Create a non-root user to use if preferred - see https://aka.ms/vscode-remote/containers/non-root-user.
# [Optional] Add sudo support for the non-root user
# [hack] Set disable_coredump false for Ubuntu Focal
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && apt-get install -y sudo lsb-release \
    && if [ "$(lsb_release -cs)" = "focal" ]; then echo "Set disable_coredump false" > /etc/sudo.conf; fi \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Clean up
RUN apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog


# Setup demo environment variables
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    LC_ALL=C.UTF-8 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

USER vscode
ENV HOME /home/vscode
WORKDIR $HOME

RUN git clone --branch ${BRANCH} https://github.com/feelpp/mqs.git; \
    mkdir /home/vscode/build && cd /home/vscode/build; \
    echo "cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-9 -DCMAKE_INSTALL_PREFIX=./install ../mqs; " \
    echo "make -j3;" \
    echo "[sudo] make install;" 

CMD ["/bin/bash"]