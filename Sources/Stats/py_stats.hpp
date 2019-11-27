#ifndef NQS_PYSTATS_HPP
#define NQS_PYSTATS_HPP

#include <pybind11/pybind11.h>

#include "common_types.hpp"
#include "obs_manager.hpp"

namespace py = pybind11;

namespace nqs {

namespace detail {

py::dict GetItem(const ObsManager& self, const std::string& name) {
  py::dict dict;
  self.InsertAllStats(name, dict);
  return dict;
}

}  // namespace detail

void AddStatsModule(py::module& m) {
  auto subm = m.def_submodule("stats");

  py::class_<ObsManager>(subm, "ObsManager")
      .def("__getitem__", &detail::GetItem, py::arg("name"))
      .def("__getattr__", &detail::GetItem, py::arg("name"))
      .def("__contains__", &ObsManager::Contains, py::arg("name"))
      .def("__len__", &ObsManager::Size)
      .def("keys", &ObsManager::Names)
      .def("__repr__", [](const ObsManager& self) {
        std::string s("<nqs.stats.ObsManager: size=");
        auto size = self.Size();
        s += std::to_string(size);
        if(size > 0) {
          s += " [";
          for (const auto& name : self.Names()) {
            s += name + ", ";
          }
          // remove last comma + space:
          s.pop_back();
          s.pop_back();

          s += "]";
        }
        return s + ">";
      });
}
}  // namespace nqs

// Expose the Stats object to Python as dict
namespace pybind11 {
namespace detail {

using NkStatsType = nqs::Binning<double>::Stats;

template <>
struct type_caster<NkStatsType> {
 public:
  PYBIND11_TYPE_CASTER(NkStatsType, _("_Stats"));

  static handle cast(NkStatsType src, return_value_policy /* policy */,
                     handle /* parent */) {
    py::dict dict;
    dict["Mean"] = src.mean;
    dict["Sigma"] = src.sigma;
    dict["Taucorr"] = src.taucorr;
    return dict.release();
  }
};

}  // namespace detail
}  // namespace pybind11

#endif  // NQS_PYSTATS_HPP
