FROM feelpp/feelpp-toolboxes:feature-biotsavart-v0.108.0-beta.1-ubuntu-20.04

USER root

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

ARG BRANCH=biotsavart-bcs

# Setup demo environment variables
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    LC_ALL=C.UTF-8 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

# install Feelpp from BinTray Debian repository
RUN apt update && \
    apt install -y lsb-release sudo && \
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
    
RUN apt install -y emacs-nox vim-nox

# Quick and Dirty hack for Focal (see https://bugzilla.redhat.com/show_bug.cgi?id=1773148)
RUN if [ $(lsb_release -cs) = "focal" ]; then sudo echo "Set disable_coredump false" > /etc/sudo.conf; fi

# #set up user so that we do not run as root
# # Changing the password does not work on certain debian and is not needed
# # echo "feelpp:docker" | chpasswd && \
# RUN useradd -m -s /bin/bash -G sudo,video feelpp && \
#      mkdir -p  /etc/sudoers.d/ && \
#      echo "feelpp ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers.d/feelpp


# get feelpp version
RUN version=$(feelpp_mesh_partitioner --version | grep version | cut -d " " -f4) && \
    echo "Feelpp: $version"

USER feelpp
# ENV HOME /home/feelpp
WORKDIR $HOME

RUN cd $HOME; pwd; ls; \
    cd src; git clone --depth 1 --branch $BRANCH https://github.com/feelpp/mqs.git; \
    cd; mkdir -p build && cd build; \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-9 $HOME/src/mqs; \
    make -j3;\
    sudo make install; \
    echo "mkdir check && cd check"; \
    echo "mpirun --bind-to core --mca btl vader,self -np 4 feelpp_mqs_form --config-file /usr/local/share/feelpp/data/testcases/mqs/cases" ; \
    make clean;

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog


CMD ["/bin/bash"]