
#include <fstream>
#include <string>
#include <vector>
#include "Utils/json_utils.hpp"

std::vector<nqs::json> GetObservableInputs() {
  std::vector<nqs::json> input_tests;
  nqs::json pars;

  std::vector<std::vector<double>> sx = {{0, 1}, {1, 0}};
  std::vector<std::vector<double>> szsz = {
      {1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, 1}};
  Complex Iu(0, 1);
  std::vector<std::vector<Complex>> sy = {{0, Iu}, {-Iu, 0}};

  pars.clear();
  pars["Hilbert"]["QuantumNumbers"] = {1, -1};
  pars["Hilbert"]["Size"] = 10;

  nqs::json parsobs;
  parsobs["Name"] = "Observable_1";
  parsobs["Operators"] = {sx, szsz, szsz, sx,   sy, sy,
                          sy, szsz, sx,   szsz, sy, szsz};
  parsobs["ActingOn"] = {{0}, {0, 1}, {1, 0}, {1},    {2}, {3},
                         {4}, {4, 5}, {5},    {6, 8}, {9}, {7, 0}};

  pars["Observables"].push_back(parsobs);

  input_tests.push_back(pars);
  return input_tests;
}
