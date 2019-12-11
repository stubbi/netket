NQS is built on top of [NetKet](https://github.com/netket/netket) for the classical simulation of Quantum Circuits using Restricted Boltzmann machines.
It is a Python library built on C++ primitives.

*Important:* The code is currently optimized for the `skylake-avx512` architecture. In order to change that, replace the compiler arguments in the [CmakeLists.txt](https://github.com/stubbi/nqs/blob/master/CMakeLists.txt#L112).