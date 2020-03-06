Bootstrap: docker
From: python:3.7

%files
    ./ nqs/

%post
    apt-get update
    apt-get install -y --no-install-recommends apt-transport-https ca-certificates gnupg software-properties-common wget
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add -
    apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
    apt-get install -y --no-install-recommends cmake openmpi-bin libopenmpi-dev libatlas-base-dev python-dev python-pip libssl-dev
    pip install -U pip setuptools numpy scipy
    pip install matplotlib
    pip install pandas
    chmod -R 755 /nqs/
    cd nqs
    python setup.py install