#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <algorithm>

const int nodes_count = 18;

double func_value(const double& x) {
  return cos(x) * x;
}

double func_deriv_value(const double& x) {
  return -sin(x) * x + cos(x);
}

std::vector<double> solve_system(std::vector<std::vector<double>>& matrix, std::vector<double>& B, int nodes_count) {
  std::vector<int> this_permutation(nodes_count);

  for(int i = 0; i < nodes_count; ++i) {
    this_permutation[i] = i;
  }

  for(int i = 0; i < nodes_count - 1; ++i) {
    //finding max-abs elem
    int max_index = i;
    double max_elem = std::abs(matrix[i][i]);
    for(int j = i + 1; j < nodes_count; ++j) {
      double potential_max = std::abs(matrix[j][i]);
      if (potential_max > max_elem) {
        max_elem = potential_max;
        max_index = j;
      }
    }

    // swapping rows
    if (i != max_index) {
      std::swap(matrix[i], matrix[max_index]);
      std::swap(this_permutation[i], this_permutation[max_index]);
    }

    //filling L[i][i] elem
    if (matrix[i][i] == 0) {
      std::cout << "TRUBA" << std::endl;
      return {};
    }

    //filling L matrix and doing Gauss
    for(int j = i + 1; j < nodes_count; ++j) {
      if (matrix[j][i] != 0) {
        double mul_value = matrix[j][i] / matrix[i][i];
        for(int t = i + 1; t < nodes_count; ++t) {
          matrix[j][t] -= mul_value * matrix[i][t];
        }
      }
    }

    for(int j = i + 1; j < nodes_count; ++j) {
      matrix[i][j] /= matrix[i][i];
    }
  }

//  DY_1=B
  for(int i = 0; i < nodes_count; ++i) {
    if (this_permutation[i] != i) {
      std::swap(B[this_permutation[this_permutation[i]]], B[this_permutation[i]]);
      std::swap(this_permutation[this_permutation[i]], this_permutation[i]);
    }
  }

  //LY_2 = B'
  for(int i = 0; i < nodes_count; ++i) {
    if (std::abs(B[i]) > 0.00001) {
      if (matrix[i][i] == 0) {
        std::cout << "Degenerate system";
        return {};
      }
      B[i] /= matrix[i][i];
      //matrix[i][i] = 1;
      for(int j = i + 1; j < nodes_count; ++j) {
        B[j] -= B[i] * matrix[j][i];
      }
    }
  }

  //UX = B''
  for(int i = nodes_count - 1; i >= 0; --i) {
    if (std::abs(B[i]) > 0.00001) {
      for(int j = i - 1; j >= 0; --j) {
        B[j] -= B[i] * matrix[j][i];
      }
    }
  }

  return B;
}

std::vector<std::pair<double, int>> task_1_2(std::vector<double>& nodes) {
  int nodes_count = nodes.size();
  clock_t start_time = clock();
  std::vector<std::vector<double>> matrix(nodes_count,
                                          std::vector<double>(nodes_count));
  std::vector<double> B(nodes_count);
  for(int i = 0; i < nodes_count; ++i) {
    for(int j = 0; j < nodes_count; ++j) {
      matrix[i][j] = std::pow(nodes[i], nodes_count - j - 1);
    }
    B[i] = func_value(nodes[i]);
  }
  B = solve_system(matrix, B, nodes_count);
  clock_t end_time = clock();

  long double search_time = (long double) (end_time - start_time) / CLOCKS_PER_SEC;
//  std::cout << "Time: " << search_time << std::endl;

  std::vector<std::pair<double, int>> ret;
  for(int i = 0; i < nodes_count; ++i) {
    ret.emplace_back(B[i], nodes_count - i - 1);
//    if (i == nodes_count - 1) {
//      std::cout << std::fixed << std::setprecision(12) << "(" << B[i] << ")*x^" << nodes_count - i - 1 << std::endl;
//    } else {
//      std::cout << std::fixed << std::setprecision(12) << "(" << B[i] << ")*x^" << nodes_count - i - 1 << "+";
//    }
  }
  return ret;
}

void task_3_4(std::vector<double>& nodes) {
  clock_t start_time = clock();
  auto temp = nodes;
  for(const auto& elem : temp) {
    nodes.push_back(elem);
  }
  std::sort(nodes.begin(), nodes.end());
  std::vector<std::vector<double>> table(nodes.size());
  int rows = nodes.size();
  int cols = rows + 1;
  for(int i = 0; i < rows; ++i)
    table[i].resize(2 + i);
  for(int i = 0; i < rows; ++i) {
    table[i][0] = nodes[i];
    table[i][1] = func_value(nodes[i]);
  }

  int sub_distance = 1;
  for(int i = 1; i < cols - 1; ++i) {
    for(int j = i - 1; j < rows - 1; ++j) {
      double value = (table[j + 1][i] - table[j][i]) /
          (table[j - i + 1 + sub_distance][0] - table[j - i + 1][0]);
      if (fabs(table[j - i + 1 + sub_distance][0] - table[j - i + 1][0]) < 1e-8) {
        value = func_deriv_value(nodes[j - 1 + 1]);
      }
      table[j + 1][i + 1] = value;
    }
    ++sub_distance;
  }

  clock_t end_time = clock();
  long double search_time = (long double) (end_time - start_time) / CLOCKS_PER_SEC;
  std::cout << "Time: " << search_time << std::endl;

  //print
  std::cout << table[0][1];
  for(int i = 2; i <= rows; ++i) {
    std::cout << std::fixed << std::setprecision(5) << "+(" << table[i - 1][i] << ")*";
    for(int j = 0; j < i - 1; ++j) {
      if (j != i - 2) {
        std::cout << std::fixed << std::setprecision(5) << "(x-" << table[j][0] << ")*";
      } else {
        std::cout << std::fixed << std::setprecision(5) << "(x-" << table[j][0] << ")";
      }
    }
  }
}

std::vector<double> get_equally_spaced_nodes(int count) {
  std::vector<double> nodes(count);
  double step = 22.0 / (count - 1);
  double start = 0;
  for(int i = 0; i < count; ++i) {
    nodes[i] = start + i * step;
  }
  return nodes;
}

std::vector<double> get_cheb_nodes(int count) {
  std::vector<double> nodes(count);
  for(int i = 0; i < count; ++i) {
    nodes[i] = 11 + 11 * cos(3.14 * (2 * i + 1) / (2 * count));
  }
  return nodes;
}

double get_polynom_value(std::vector<std::pair<double, int>>& coeff_power, double x) {
  double res = 0;
  for(const auto& elem : coeff_power) {
    res += elem.first * std::pow(x, elem.second);
  }
  return res;
}

double task_5() {
  const double value = std::sqrt(2) / 2;
//  const double value = M_PI;
  auto nodes = get_equally_spaced_nodes(100);
  int nearest_ind = 0;
  double nearest_value = INT8_MAX;
  for(int i = 0; i < 100; ++i) {
    if (fabs(nodes[i] - value) < nearest_value) {
      nearest_value = fabs(nodes[i] - value);
      nearest_ind = i;
    }
  }

  std::vector<double> nodes_to_check = {nodes[nearest_ind]};
  std::vector<std::pair<double, int>> pol;
  double cur_value = 1;
  double prev_value = 0;
  bool right_margin = false;
  bool left_margin = false;
  int left_ind = nearest_ind;
  int right_ind = nearest_ind;
  int nodes_count = 1;
  while ((!left_margin && !right_margin) || fabs(cur_value - prev_value) >= 1e-8) {
    prev_value = cur_value;
    if (right_ind > 100) right_margin = true;
    if (left_ind < 0) left_margin = true;
    if (!left_margin) {
      --left_ind;
      nodes_to_check.push_back(nodes[left_ind]);
      ++nodes_count;
    }
    if (!right_margin) {
      ++right_ind;
      nodes_to_check.push_back(nodes[right_ind]);
      ++nodes_count;
    }
    pol = task_1_2(nodes_to_check);
    cur_value = get_polynom_value(pol, value);

  }
  std::cout << "Nodes count: " << nodes_count << "Value: " <<  cur_value << std::endl;
  return cur_value;
}

void task_6(std::vector<double>& nodes) {
  int N = nodes_count - 1;
  std::vector<double>
      h(N + 1),
      x(N + 1),
      l(N + 1),
      b(N + 1),
      y(N + 1),
      delta(N + 1),
      lambda(N + 1),
      c(N + 1),
      d(N + 1);
  for(int i = 0; i <= N; ++i) {
    x[i] = nodes[i];
    y[i] = func_value(nodes[i]);
  }

  for(int k = 1; k <= N; k++) {
    h[k] = x[k] - x[k - 1];
    if (h[k] == 0) {
      std::cout << "Error: x[" << k << "]=x[" << k - 1 << ']';
      return;
    }
    l[k] = (y[k] - y[k - 1]) / h[k];
  }
  delta[1] = -h[2] / (2 * (h[1] + h[2]));
  lambda[1] = 1.5 * (l[2] - l[1]) / (h[1] + h[2]);
  for(int k = 3; k <= N; k++) {
    delta[k - 1] = -h[k] / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]);
    lambda[k - 1] = (3 * l[k] - 3 * l[k - 1] - h[k - 1] * lambda[k - 2]) /
        (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]);
  }
  c[0] = 0;
  c[N] = 0;
  for(int k = N; k >= 2; k--) {
    c[k - 1] = delta[k - 1] * c[k] + lambda[k - 1];
  }
  for(int k = 1; k <= N; k++) {
    d[k] = (c[k] - c[k - 1]) / (3 * h[k]);
    b[k] = l[k] + (2 * c[k] * h[k] + h[k] * c[k - 1]) / 3;
  }
  //printresult(N, y, b, c, d);
  for(int i = 1; i <= N; ++i) {
    std::cout << "F_" << i << " = " << y[i] <<
              "+(" << b[i] << ")*(x-(" << nodes[i] << "))+(" <<
              c[i] << ")*(x-(" << nodes[i] << "))^2+(" <<
              d[i] << ")*(x-(" << nodes[i] << "))^3" << std::endl;
  }
}

double scalar_mul(std::vector<double>& l, std::vector<double>& r) {
  double ret = 0;
  for(int i = 0; i < l.size(); ++i) {
    ret += l[i] * r[i];
  }
  return ret;
}

void task_7(int n) {
  std::vector<double> nodes;
  srand(time(nullptr));
  for(int i = 0; i < 100; ++i) {
    double value = (rand() % (22001)) / 1000.0;
    nodes.push_back(value);
  }
  std::sort(nodes.begin(), nodes.end());

  std::vector<std::vector<double>> table(n + 1, std::vector<double>(100));
  std::vector<double> f(100);
  for(int i = 0; i < 5; ++i) {
    for(int j = 0; j < n + 1; ++j) {
      table[j][i] = std::pow(nodes[i], j);
    }
    f[i] = func_value(nodes[i]);
  }

  std::vector<std::vector<double>> matrix(n + 1, std::vector<double>(n + 1));
  std::vector<double> b(n + 1);

  for(int i = 0; i < n + 1; ++i) {
    for(int j = 0; j < n + 1; ++j) {
      matrix[i][j] = scalar_mul(table[i], table[j]);
    }
    b[i] = scalar_mul(f, table[i]);
  }

  auto ans = solve_system(matrix, b, n + 1);

  for(int i = 0; i < n + 1; ++i) {
    if (i != n) {
      std::cout << std::fixed << std::setprecision(8) << '(' << ans[i] << ")*x^" << i << '+';
    } else {
      std::cout << std::fixed << std::setprecision(8) << '(' << ans[i] << ")*x^" << i;
    }
  }
}

int main() {
//  /*Пункт 1,3*/auto nodes = get_equally_spaced_nodes(nodes_count);
   /*Пункт 2,4*/auto nodes = get_cheb_nodes(nodes_count);
//  task_1_2(nodes);
//  task_3_4(nodes);
  task_5();
//  task_6(nodes);
//  task_7(6);

  return 0;
}