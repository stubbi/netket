NQS is built on top of [NetKet](https://github.com/netket/netket) for the classical simulation of Quantum Circuits using Restricted Boltzmann machines.
It is a Python library built on C++ primitives.

# Usage

The code is written in C++ and utilizes MPI to allow fast execution. A Python interface is provided through PyBind. You can either clone this repository and install the sources on your own using:

```
    git clone https://github.com/stubbi/nqs.git
    cd nqs
    python setup.py install
```

*Important:* The code is currently optimized for the `skylake-avx512` architecture. In order to change that, replace the compiler arguments in the [CmakeLists.txt](https://github.com/stubbi/nqs/blob/master/CMakeLists.txt#L112).

There is also the alternative to use the [Singularity Container](https://singularity-hub.org/) which has the package preinstalled and run your scripts inside them. For this, you need to have singularity installed on your machine:

```
    singularity pull shub://stubbi/nqs
```

Afterwards, you can simply add `nqs` as a dependency to your python project:

```
import nqs
```
