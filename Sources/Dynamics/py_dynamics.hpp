#ifndef NQS_DYNAMICS_PY_DYNAMICS_HPP
#define NQS_DYNAMICS_PY_DYNAMICS_HPP

#include <pybind11/pybind11.h>

namespace nqs {

void AddDynamicsModule(pybind11::module m);

}  // namespace nqs

#endif  // NQS_DYNAMICS_PY_DYNAMICS_HPP
