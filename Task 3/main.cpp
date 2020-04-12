#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

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
  for(int i = 0; i < nodes_count; ++i) {
    for(int j = 0; j < nodes_count; ++j) {
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
      for(int k = 0; k < ans.size(); ++k) {
        ans[k] += x_pol[k / nodes_count] * y_pol[k % nodes_count];
      }
    }

  for(int i = 0; i < nodes_count; ++i)
    for(int j = 0; j < nodes_count; ++j) {
      std::cout << std::fixed << std::setprecision(3);
      if (fabs(ans[nodes_count * i + j]) <= 1e-3) continue;
      if (nodes_count * i + j != ans.size() - 1) {
        std::cout << ans[nodes_count * i + j];
        if (i >= 10) std::cout << "*x^(" << i << ")";
        else std::cout << "*x^" << i << "";
        if (j >= 10) std::cout << "y^(" << j << ")";
        else std::cout << "y^" << j;
        if (nodes_count * i + j != ans.size() - 1 && ans[nodes_count * i + j + 1] >= 1e-4) {
          std::cout << '+';
        }
      }
    }
}

std::vector<double> vector_sum(const std::vector<double>& lhs,
                               const std::vector<double>& rhs) {
  std::vector<double> res(lhs.size());
  for(int i = 0; i < lhs.size(); ++i) {
    res[i] = lhs[i] + rhs[i];
  }
  return res;
}

std::vector<double> three_diagonal(std::vector<double>& matrix,
                                   std::vector<double>& b) {
  for(int i = 0; i < matrix.size() - 5; i += 3) {
    if (fabs(matrix[i]) < fabs(matrix[i + 3])) {
      std::swap(matrix[i], matrix[i + 3]);
      std::swap(matrix[i + 1], matrix[i + 4]);
      std::swap(matrix[i + 2], matrix[i + 5]);
      std::swap(b[(i / 3)], b[(i / 3) + 1]);
    }
    double k = -1 * (matrix[i + 3] / matrix[i]);
    matrix[i + 3] = matrix[i + 4] + k * matrix[i + 1];
    matrix[i + 4] = matrix[i + 5] + k * matrix[i + 2];
    matrix[i + 5] = 0;
    b[(i / 3) + 1] += k * b[i / 3];
  }
  b[b.size() - 1] /= matrix[matrix.size() - 3];
  b[b.size() - 2] -= matrix[matrix.size() - 5] * b[b.size() - 1];
  b[b.size() - 2] /= matrix[matrix.size() - 6];
  for(int i = b.size() - 3; i >= 0; i--) {
    b[i] += b[i + 1] * (-1) * matrix[3 * i + 1];
    b[i] += b[i + 2] * (-1) * matrix[3 * i + 2];
    b[i] /= matrix[3 * i];
  }
  return b;
}

void three_diag_system(const double& val,
                       std::vector<double>& u,
                       std::vector<std::vector<double>>& M,
                       int j) {
  std::vector<double> matrix;
  double a = val / 6;
  double b = -2 * val / 3;
  double c = val / 6;
  matrix.push_back(-b);
  matrix.push_back(c);
  matrix.push_back(0);
  for(int i = 1; i < u.size() - 1; ++i) {
    matrix.push_back(a);
    matrix.push_back(-b);
    matrix.push_back(c);
  }
  if (u.size() == 1) {
    M[j].push_back(0);
    M[j].push_back(u[0] / (-b));
    M[j].push_back(0);
  }
  if (u.size() > 1) {
    matrix.push_back(a);
    matrix.push_back(-b);
    matrix.push_back(0);
    M[j].push_back(0);
    for(auto& i : three_diagonal(matrix, u)) {
      M[j].push_back(i);
    }
    M[j].push_back(0);
  }
}

std::vector<std::vector<std::vector<double>>> task_2(std::vector<double>& x,
                                                     std::vector<double>& y) {
  double h = x[1] - x[0];
  double t = y[1] - y[0];
  std::vector<double> ans(16);
  std::vector<std::vector<std::vector<double>>>
      res((nodes_count - 1), std::vector<std::vector<double>>(nodes_count - 1));
  std::vector<std::vector<std::vector<double>>>
      deltas(nodes_count, std::vector<std::vector<double>>(nodes_count)),
      polynoms(nodes_count, std::vector<std::vector<double>>(nodes_count));

  std::vector<std::vector<double>> M(y.size());
  for(int j = 0; j < y.size(); j++) {
    std::vector<double> d(x.size() - 2);
    for(int i = 1; i < x.size() - 1; i++) {
      d[i - 1] = (func_value(x[i + 1], y[j]) - func_value(x[i], y[j])) / h -
          (func_value(x[i], y[j]) - func_value(x[i - 1], y[j])) / h;
    }
    three_diag_system(h, d, M, j);
  }

  std::vector<std::vector<double>> K(x.size());
  for(int i = 0; i < x.size(); ++i) {
    std::vector<double> d(y.size() - 2);
    for(int j = 1; j < y.size() - 1; ++j) {
      d[j - 1] = (M[j + 1][i] - M[j][i]) / t - (M[j][i] - M[j - 1][i]) / t;
    }
    three_diag_system(t, d, K, i);
  }

  std::vector<std::vector<double>> L(x.size());
  for(int i = 0; i < x.size(); ++i) {
    std::vector<double> d(y.size() - 2);
    for(int j = 1; j < y.size() - 1; ++j) {
      d[j - 1] = (func_value(x[i], y[j + 1]) - func_value(x[i], y[j])) / t
          - (func_value(x[i], y[j]) - func_value(x[i], y[j - 1])) / t;
    }
    three_diag_system(t, d, L, i);
  }

  for(int i = 0; i < x.size(); ++i) {
    for(int j = 1; j < y.size(); ++j) {
      std::vector<double> splayin_y = {0};
      std::vector<double> delta = {0};
      std::vector<double> polynom = {y[j], -1};
      std::vector<double> s = polynom;
      std::vector<double> c = polynom;
      polynom = mul_pols(polynom, polynom);
      polynom = mul_pols(polynom, s);
      std::vector<double> d = polynom;
      mul_pol_by_scalar(d, L[i][j - 1] / (6 * t));
      delta = vector_sum(delta, d);
      mul_pol_by_scalar(polynom, K[i][j - 1] / (6 * t));
      splayin_y = vector_sum(splayin_y, polynom);
      mul_pol_by_scalar(c, func_value(x[i], y[j - 1]) / t - L[i][j - 1] * t / 6);
      mul_pol_by_scalar(s, M[j - 1][i] / t - L[i][j - 1] * t / 6);
      splayin_y = vector_sum(splayin_y, s);
      delta = vector_sum(delta, c);
      polynom = {0};
      polynom = {-y[j - 1], 1};
      s = polynom;
      c = polynom;
      polynom = mul_pols(polynom, polynom);
      polynom = mul_pols(polynom, s);
      d = polynom;
      mul_pol_by_scalar(d, L[i][j] / (6 * t));
      mul_pol_by_scalar(polynom, K[i][j] / (6 * t));
      splayin_y = vector_sum(splayin_y, polynom);
      delta = vector_sum(delta, d);
      mul_pol_by_scalar(s, M[j][i] / t - L[i][j] * t / 6);
      mul_pol_by_scalar(c, func_value(x[i], y[j]) / t - L[i][j] * t / 6);
      splayin_y = vector_sum(splayin_y, s);
      delta = vector_sum(delta, c);

      deltas[i][j] = delta;
      polynoms[i][j] = splayin_y;
    }
  }

  for(int i = 1; i < x.size(); ++i) {
    for(int j = 1; j < y.size(); ++j) {
      for(auto& r : ans) {
        r = 0;
      }
      std::vector<double> polynom_x = {x[i], -1};
      std::vector<double> s = polynom_x;
      polynom_x = mul_pols(polynom_x, s);
      polynom_x = mul_pols(polynom_x, s);
      s.push_back(0);
      s.push_back(0);
      mul_pol_by_scalar(polynom_x, 1 / (6 * h));
      for(int k = 0; k < ans.size(); ++k) {
        ans[k] += polynom_x[k / 4] * polynoms[i - 1][j][k % 4];
      }

      polynom_x = {-x[i - 1], 1};
      std::vector<double> c = polynom_x;
      polynom_x = mul_pols(polynom_x, c);
      polynom_x = mul_pols(polynom_x, c);
      c.push_back(0);
      c.push_back(0);
      mul_pol_by_scalar(polynom_x, 1 / (6 * h));
      for(int k = 0; k < ans.size(); ++k) {
        ans[k] += polynom_x[k / 4] * polynoms[i][j][k % 4];
      }

      std::vector<double> polynom_y = deltas[i][j];
      std::vector<double> vector_suming_y = polynoms[i][j];
      mul_pol_by_scalar(vector_suming_y, -pow(h, 2) / 6);
      polynom_y = vector_sum(polynom_y, vector_suming_y);
      mul_pol_by_scalar(c, 1 / h);
      for(int k = 0; k < ans.size(); ++k) {
        ans[k] += c[k / 4] * polynom_y[k % 4];
      }

      polynom_y = deltas[i - 1][j];
      vector_suming_y = polynoms[i - 1][j];
      mul_pol_by_scalar(vector_suming_y, -pow(h, 2) / 6);
      polynom_y = vector_sum(polynom_y, vector_suming_y);
      mul_pol_by_scalar(s, 1 / h);
      for(int k = 0; k < ans.size(); ++k) {
        ans[k] += s[k / 4] * polynom_y[k % 4];
      }
      res[i - 1][j - 1] = ans;
    }
  }
  return res;
}

int main() {
  std::ofstream out("output1.txt");
  auto x = get_equally_spaced_nodes();
  auto y = get_equally_spaced_nodes();
//  task_1(x, y);
  auto splain = task_2(x, y);
  for(int i = 0; i < x.size() - 1; ++i) {
    for(int j = 0; j < y.size() - 1; ++j) {
      auto s = splain[i][j];
//      out << "Splain with x in [" << x[i] << ", " << x[i + 1] << "]; y in [" << y[j] << ", " << y[j + 1] << "]" << std::endl;
      out << "def makePolynom_" << i << "_" << j << "() :" << std::endl;
      out << "  a = np.arange ( " << x[i] << ", " << x[i + 1] << ", 0.1)" << std::endl;
      out << "  b = np.arange ( " << y[j] << ", " << y[j + 1] << ", 0.1)" << std::endl;
      out << "  x, y = np.meshgrid(a, b)" << std::endl;
      out << "  z = ";
      for(int m = 0; m < 3; m++) {
        for(int k = 0; k < 4; ++k) {
          out << s[k] << " * ( x ** " << m << ") * ( y ** " << k << ") + ";
        }
      }
      for(int k = 0; k < 3; ++k) {
        out << s[k] << " * ( x ** " << 4 << ") * ( y ** " << k << ") + ";
      }
      out << s[4] << " * ( x ** " << 4 << ") * ( y ** " << 4 << ")";
      out << std::endl;
      out << "  return x, y, z" << std::endl;
      out << std::endl;
    }
  }
  return 0;
}
