ARG BASE_CONTAINER=ubuntu:focal
FROM $BASE_CONTAINER

LABEL maintainer="Christophe Trophime <christophr.trophime@lncmi.cnrs.fr>"

USER root

# select feelpp version: stable or latest
ARG FEELPP_CHANNEL=stable

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Setup demo environment variables
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    LC_ALL=C.UTF-8 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

# use proposed to get missing packages (not really SAFE)
RUN apt update && \
    apt install -y lsb-release wget sudo && \
    echo "deb http://archive.ubuntu.com/ubuntu $(lsb_release -cs)-proposed main restricted"  |  tee -a /etc/apt/sources.list && \
    echo "deb http://archive.ubuntu.com/ubuntu $(lsb_release -cs)-proposed universe"  |  tee -a /etc/apt/sources.list && \
    echo "deb http://archive.ubuntu.com/ubuntu $(lsb_release -cs)-proposed multiverse"  |  tee -a /etc/apt/sources.list && \
    apt update

# install Feelpp from BinTray Debian repository
RUN apt update && \
    apt install -y lsb-release sudo wget gnupg && \
    echo "deb https://dl.bintray.com/feelpp/ubuntu $(lsb_release -cs) $FEELPP_CHANNEL" > /etc/apt/sources.list.d/feelpp.list && \
    wget -qO - https://bintray.com/user/downloadSubjectPublicKey?username=bintray | apt-key add - && \
    apt update && \
    apt-cache search feelpp && \
    echo "FEELPP_CHANNEL=${FEELPP_CHANNEL}" && \
    if [ "$FEELPP_CHANNEL" = "stable" ]; then \
       apt install -y libfeelpp-toolboxes-dev feelpp-tools; \
    elif [ "$FEELPP_CHANNEL" = "latest" ]; then \
       apt install -y libfeelpp-toolboxes1-core-dev feelpp-tools; \
     else { \
       echo "FEELPP_CHANNEL unknown: ${FEELPP_CHANNEL}"; \
       exit 1; \
    } fi && \
    apt install -y clang-9 g++ git cmake

RUN apt install -y emacs-nox vim-nox

# Quick and Dirty hack for Focal (see https://bugzilla.redhat.com/show_bug.cgi?id=1773148)
RUN if [ $(lsb_release -cs) = "focal" ]; then sudo echo "Set disable_coredump false" > /etc/sudo.conf; fi

#set up user so that we do not run as root
# Changing the password does not work on certain debian and is not needed
# echo "feelpp:docker" | chpasswd && \
RUN useradd -m -s /bin/bash -G sudo,video feelpp && \
     mkdir -p  /etc/sudoers.d/ && \
     echo "feelpp ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers.d/feelpp

# get feelpp version
RUN version=$(feelpp_toolbox_thermoelectric --version | grep version | cut -d " " -f4) && \
    echo "Feelpp: $version"

USER feelpp
ENV HOME /home/feelpp
WORKDIR $HOME

RUN cd $HOME; pwd; git clone https://github.com/feelpp/mqs.git; \
    mkdir build && cd build; \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-9 ../mqs; \
    make -j3;\
    sudo make install; \
    echo "mkdir check && cd check"; \
    echo "mpirun --bind-to core --mca btl vader,self -np 4 feelpp_mqs_form --config-file /usr/local/share/feelpp/data/testcases/mqs/cases" ; \
    make clean;

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog



