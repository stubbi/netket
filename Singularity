Bootstrap: docker
From: ubuntu:18.04

%files
    ./

%post
    apt-get update
    apt-get install -y --no-install-recommends apt-transport-https ca-certificates gnupg software-properties-common wget
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add -
    apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
    apt-get install -y --no-install-recommends cmake libmpich-dev mpich libatlas-base-dev python-pip libssl-dev
    python setup.py install

