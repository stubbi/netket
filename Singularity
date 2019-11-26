Bootstrap: docker
From: ubuntu:18.04

%files
    ./*

%post
    apt-get update
    apt-get install -y --no-install-recommends cmake libmpich-dev mpich libatlas-base-dev python-pip libssl-dev
    python setup.py install

