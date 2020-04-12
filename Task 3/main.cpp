#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

const int nodes_count = 18;

double func_value(const double& x, const double& y) {
  return cos(x * y);
}

std::vector<double> get_equally_spaced_nodes() {
  std::vector<double> nodes(nodes_count);
  double step = 4.0 / (nodes_count - 1);
  double start = 0;
  for(int i = 0; i < nodes_count; ++i) {
    nodes[i] = start + i * step;
  }
  return nodes;
}

std::vector<double> mul_pols(const std::vector<double>& lhs, const std::vector<double>& rhs) {
  std::vector<double> res((lhs.size() - 1) + rhs.size());
  for(int i = 0; i < lhs.size(); ++i)
    for(int j = 0; j < rhs.size(); ++j)
      res[i + j] += rhs[j] * lhs[i];
  return res;
}

void mul_pol_by_scalar(std::vector<double>& pol, const double& scalar) {
  for(auto& elem : pol)
    elem *= scalar;
}

double pol_value(const std::vector<double>& pol, const double& x) {
  double sum = 0;
  for(int i = 0; i < pol.size(); ++i)
    sum += pow(x, i) * pol[i];
  return sum;
}

void task_1(std::vector<double>& x, std::vector<double>& y) {
  std::vector<double> ans(nodes_count * nodes_count);
  std::vector<std::vector<std::vector<double>>>
      w_x(nodes_count, std::vector<std::vector<double>>(nodes_count)),
      w_y(nodes_count, std::vector<std::vector<double>>(nodes_count));
  for(int i = 0; i < nodes_count; i++) {
    for(int j = 0; j < nodes_count; j++) {
      w_x[i][j] = {-x[j], 1};
      w_y[i][j] = {-y[j], 1};
    }
    w_x[i][i] = {1};
    w_y[i][i] = {1};
  }

  for(int i = 0; i < nodes_count; ++i)
    for(int j = 0; j < nodes_count; ++j) {
      std::vector<double> x_pol = {1};
      std::vector<double> y_pol = {1};
      for(auto& elem : w_x[i])
        x_pol = mul_pols(x_pol, elem);
      for(auto& elem : w_y[j])
        y_pol = mul_pols(y_pol, elem);

      double z = pol_value(x_pol, x[i]);
      double q = pol_value(y_pol, y[j]);
      mul_pol_by_scalar(x_pol, func_value(x[i], y[j]) / (z * q));
      for(int k = 0; k < ans.size(); k++) {
        ans[k] += x_pol[k / nodes_count] * y_pol[k % nodes_count];
      }
    }

  for(int i = 0; i < nodes_count; ++i)
    for(int j = 0; j < nodes_count; ++j) {
      std::cout << std::fixed << std::setprecision(5);
      if (fabs(ans[nodes_count * i + j]) <= 1e-5) continue;
      if (nodes_count * i + j != ans.size() - 1) {
        std::cout << '(' << ans[nodes_count * i + j] << ")*x^(" << i << ")y^(" << j << ")+";
      } else {
        std::cout << '(' << ans[nodes_count * i + j] << ")*x^(" << i << ")y^(" << j << ")";
      }
    }
}

int main() {
  auto x = get_equally_spaced_nodes();
  auto y = get_equally_spaced_nodes();
  task_1(x, y);
  std::vector<double> a = {1, 4, 2};
  std::vector<double> b = {0, 2, 7};
  return 0;
}
