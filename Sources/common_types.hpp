#ifndef NQS_COMMON_TYPES_HPP
#define NQS_COMMON_TYPES_HPP
/**
 * This header contains standard type aliases to be used throughout the Nqs
 * codebase.
 */

#include <complex>
#include <cstddef>

#include <Eigen/Core>

namespace nqs {

using Index = std::ptrdiff_t;
using Complex = std::complex<double>;

using VectorXd = Eigen::VectorXd;
using VectorXcd = Eigen::VectorXcd;
using MatrixXd = Eigen::MatrixXd;
using MatrixXcd = Eigen::MatrixXcd;

}  // namespace nqs

#endif  // NQS_COMMON_TYPES_HPP
