#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

const double eps = 1e-6;
const double deriv = 1e-8;

double norm(std::pair<double, double>& x_1, std::pair<double, double>& x_2) {
  return std::max(std::fabs(x_1.first - x_2.first),
                  std::fabs(x_1.second - x_2.second));
}

double calculate_F_1(const double& x, const double& y) {
  return x * x * y - y * y + 4 * x - 10;
}

double calculate_F_2(const double& x, const double& y) {
  return 1 / ((3 * y - 2) * (3 * y - 2)) - x - 5;
}

std::vector<double> calculate_F(const double& x, const double& y) {
  return {
      calculate_F_1(x, y),
      calculate_F_2(x, y),
  };
}

std::vector<double> calculate_minus_F(const double& x, const double& y) {
  return {
      -(calculate_F_1(x, y)),
      -(calculate_F_2(x, y)),
  };
}

std::vector<std::vector<double>> calculate_W(const double& x, const double& y) {
  return {
      {2 * x * y + 4, -2 * y + x * x},
      {-1, -6 / ((3 * y - 2) * (3 * y - 2) * (3 * y - 2))}
  };
}

double F_1_x_value(const double& x, const double& y) {
  return (calculate_F_1(x, y) - calculate_F_1(x - deriv, y)) / deriv;
}

double F_1_y_value(const double& x, const double& y) {
  return (calculate_F_1(x, y) - calculate_F_1(x, y - deriv)) / deriv;
}

double F_2_x_value(const double& x, const double& y) {
  return (calculate_F_2(x, y) - calculate_F_2(x - deriv, y)) / deriv;
}

double F_2_y_value(const double& x, const double& y) {
  return (calculate_F_2(x, y) - calculate_F_2(x, y - deriv)) / deriv;
}

int main() {

  std::pair<double, double> x_k = {INT8_MIN, INT8_MIN};
//  std::pair<double, double> x_kp = {1.7, 0.75}; // 1.92352 0.793349
//  std::pair<double, double> x_kp = {2, 0.5}; // 2.0209 0.540866
  std::pair<double, double> x_kp = {-5, 1}; // -4.76876 1.35984

  int iters = 0;
  while (norm(x_k, x_kp) >= eps) {
//    auto W = calculate_W(x_kp.first, x_kp.second);
    std::vector<std::vector<double>> W = {
        {F_1_x_value(x_kp.first, x_kp.second), F_1_y_value(x_kp.first, x_kp.second)},
        {F_2_x_value(x_kp.first, x_kp.second), F_2_y_value(x_kp.first, x_kp.second)}
    };
    auto minus_F = calculate_minus_F(x_kp.first, x_kp.second);

    W[0][1] /= W[0][0];
    minus_F[0] /= W[0][0];

    W[1][1] -= W[0][1] * W[1][0];
    minus_F[1] -= minus_F[0] * W[1][0];

    minus_F[1] /= W[1][1];
    minus_F[0] -= minus_F[1] * W[0][1];

    x_k = x_kp;
    x_kp = {x_kp.first + minus_F[0], x_kp.second + minus_F[1]};
    ++iters;
  }

  auto F = calculate_F(x_kp.first, x_kp.second);
  std::cout << iters << " iterations" << std::endl;
  std::cout << std::fixed << std::setprecision(10)
            << x_kp.first << " " << x_kp.second << std::endl;
  std::cout << std::fixed << std::setprecision(15) << F[0] << " " << F[1];

  return 0;
}