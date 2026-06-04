#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double XI = PI / 4.0;
constexpr double K1 = 1.0;
constexpr double K2 = 0.5;
constexpr double Q1 = 1.0;
constexpr double Q2 = (PI / 4.0) * (PI / 4.0);
constexpr double F1 = 1.0;
constexpr double F2 = 0.70710678118654752440;
constexpr double MU1 = 1.0;
constexpr double MU2 = 0.0;
constexpr double BETA = 1.0;

struct AnalyticConstants {
    double A{};
    double B{};
    double C{};
    double D{};
};

struct Row {
    int index{};
    double x{};
    double u{};
    double v{};
    double diff{};
};

std::string json_escape(const std::string& input) {
    std::string out;
    out.reserve(input.size() + 8);
    for (char c : input) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default: out += c; break;
        }
    }
    return out;
}

[[noreturn]] void fail(const std::string& message) {
    std::cout << "{\"error\":\"" << json_escape(message) << "\"}\n";
    std::exit(1);
}

double parse_double(const char* value, const std::string& name) {
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || *end != '\0' || !std::isfinite(parsed)) {
        fail("Invalid " + name + ": " + value);
    }
    return parsed;
}

int parse_int(const char* value, const std::string& name) {
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value || *end != '\0' || parsed < 2 || parsed > 20000) {
        fail("Invalid " + name + ": expected integer in [2, 20000]");
    }
    return static_cast<int>(parsed);
}

std::vector<double> solve_linear_4x4(double matrix[4][5]) {
    constexpr int N = 4;
    for (int col = 0; col < N; ++col) {
        int pivot = col;
        for (int row = col + 1; row < N; ++row) {
            if (std::abs(matrix[row][col]) > std::abs(matrix[pivot][col])) {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot][col]) < 1e-14) {
            throw std::runtime_error("Analytic constants system is singular");
        }
        if (pivot != col) {
            for (int j = col; j <= N; ++j) {
                std::swap(matrix[pivot][j], matrix[col][j]);
            }
        }
        const double divisor = matrix[col][col];
        for (int j = col; j <= N; ++j) {
            matrix[col][j] /= divisor;
        }
        for (int row = 0; row < N; ++row) {
            if (row == col) continue;
            const double factor = matrix[row][col];
            for (int j = col; j <= N; ++j) {
                matrix[row][j] -= factor * matrix[col][j];
            }
        }
    }

    std::vector<double> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = matrix[i][N];
    }
    return result;
}

AnalyticConstants analytic_constants() {
    const double omega = std::sqrt(2.0 * Q2);
    const double sxi = std::sinh(XI);
    const double cxi = std::cosh(XI);
    const double right_len = 1.0 - XI;
    const double sr = std::sinh(omega * right_len);
    const double cr = std::cosh(omega * right_len);
    const double particular2 = F2 / Q2;

    double m[4][5] = {
        {1.0, 0.0, 0.0, 0.0, 0.0},
        {std::cosh(XI), sxi, -1.0, 0.0, particular2 - 1.0},
        {K1 * std::sinh(XI), K1 * cxi, 0.0, -K2 * omega, 0.0},
        {0.0, 0.0, -K2 * omega * sr + BETA * cr,
         -K2 * omega * cr + BETA * sr, -BETA * particular2 + MU2}
    };

    const auto solution = solve_linear_4x4(m);
    return {solution[0], solution[1], solution[2], solution[3]};
}

double analytic_value(double x, const AnalyticConstants& constants) {
    if (x <= XI) {
        return constants.A * std::cosh(x) + constants.B * std::sinh(x) + 1.0;
    }
    const double omega = std::sqrt(2.0 * Q2);
    const double z = omega * (x - XI);
    return constants.C * std::cosh(z) + constants.D * std::sinh(z) + F2 / Q2;
}

// Harmonic average on an interval, exact for a piecewise-constant k(x).
double diffusion_on_interval(double left, double right) {
    const double length = right - left;
    if (right <= XI) return K1;
    if (left >= XI) return K2;
    return length / ((XI - left) / K1 + (right - XI) / K2);
}

// Cell averages for q and f. theta/gamma select the value at the discontinuity if it
// lies exactly on the integration cell boundary; otherwise exact piecewise integration is used.
double average_piecewise(double left, double right, double value_left, double value_right, double weight_at_break) {
    if (right <= XI) return value_left;
    if (left >= XI) return value_right;
    const double length = right - left;
    if (length <= 0.0) return weight_at_break * value_left + (1.0 - weight_at_break) * value_right;
    return ((XI - left) * value_left + (right - XI) * value_right) / length;
}

double clamp_unit(double value, const std::string& name) {
    if (value < 0.0 || value > 1.0) {
        fail(name + " must be in [0, 1]");
    }
    return value;
}

double require_positive(double value, const std::string& name) {
    if (value <= 0.0) {
        fail(name + " must be positive");
    }
    return value;
}

std::vector<double> solve_numerically(int n, double theta, double gamma) {
    const double h = 1.0 / static_cast<double>(n);
    std::vector<double> lower(n + 1, 0.0), diag(n + 1, 0.0), upper(n + 1, 0.0), rhs(n + 1, 0.0);

    diag[0] = 1.0;
    rhs[0] = MU1;

    for (int i = 1; i < n; ++i) {
        const double x = i * h;
        const double cell_left = std::max(0.0, x - 0.5 * h);
        const double cell_right = std::min(1.0, x + 0.5 * h);
        const double a_left = diffusion_on_interval((i - 1) * h, i * h);
        const double a_right = diffusion_on_interval(i * h, (i + 1) * h);
        const double q = average_piecewise(cell_left, cell_right, Q1, Q2, theta);
        const double f = average_piecewise(cell_left, cell_right, F1, F2, gamma);

        lower[i] = -a_left;
        diag[i] = a_left + a_right + q * h * h;
        upper[i] = -a_right;
        rhs[i] = f * h * h;
    }

    const double a_n = diffusion_on_interval((n - 1) * h, 1.0);
    const double q_n = Q2;
    const double f_n = F2;
    lower[n] = a_n / h;
    diag[n] = -a_n / h + BETA - 0.5 * h * q_n;
    rhs[n] = MU2 - 0.5 * h * f_n;

    std::vector<double> alpha(n + 1, 0.0), beta(n + 1, 0.0), values(n + 1, 0.0);
    if (std::abs(diag[0]) < 1e-14) throw std::runtime_error("Zero pivot at left boundary");
    alpha[0] = -upper[0] / diag[0];
    beta[0] = rhs[0] / diag[0];

    for (int i = 1; i <= n; ++i) {
        const double denom = diag[i] + lower[i] * alpha[i - 1];
        if (std::abs(denom) < 1e-14) {
            throw std::runtime_error("Zero pivot in tridiagonal solver");
        }
        alpha[i] = (i == n) ? 0.0 : -upper[i] / denom;
        beta[i] = (rhs[i] - lower[i] * beta[i - 1]) / denom;
    }

    values[n] = beta[n];
    for (int i = n - 1; i >= 0; --i) {
        values[i] = alpha[i] * values[i + 1] + beta[i];
    }
    return values;
}

void write_json_number(std::ostream& os, double value) {
    if (std::isfinite(value)) {
        os << std::setprecision(17) << value;
    } else {
        os << "null";
    }
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        fail("Usage: solver <n> <theta> <gamma> <epsilon_target>");
    }

    try {
        const int n = parse_int(argv[1], "n");
        const double theta = clamp_unit(parse_double(argv[2], "theta"), "theta");
        const double gamma = clamp_unit(parse_double(argv[3], "gamma"), "gamma");
        const double epsilon_target = require_positive(parse_double(argv[4], "epsilon_target"), "epsilon_target");
        const double h = 1.0 / static_cast<double>(n);
        const auto constants = analytic_constants();
        const auto numerical = solve_numerically(n, theta, gamma);

        std::vector<Row> rows;
        rows.reserve(n + 1);
        double epsilon1 = -1.0;
        double max_x = 0.0;
        for (int i = 0; i <= n; ++i) {
            const double x = i * h;
            const double u = analytic_value(x, constants);
            const double diff = u - numerical[i];
            const double abs_diff = std::abs(diff);
            if (abs_diff > epsilon1) {
                epsilon1 = abs_diff;
                max_x = x;
            }
            rows.push_back({i, x, u, numerical[i], diff});
        }

        std::ostringstream out;
        out << "{\"n\":" << n
            << ",\"theta\":"; write_json_number(out, theta);
        out << ",\"gamma\":"; write_json_number(out, gamma);
        out << ",\"epsilon_target\":"; write_json_number(out, epsilon_target);
        out << ",\"epsilon1\":"; write_json_number(out, epsilon1);
        out << ",\"max_x\":"; write_json_number(out, max_x);
        out << ",\"constants\":{\"A\":"; write_json_number(out, constants.A);
        out << ",\"B\":"; write_json_number(out, constants.B);
        out << ",\"C\":"; write_json_number(out, constants.C);
        out << ",\"D\":"; write_json_number(out, constants.D);
        out << "},\"points\":[";
        for (std::size_t i = 0; i < rows.size(); ++i) {
            if (i) out << ',';
            out << "{\"i\":" << rows[i].index << ",\"x\":";
            write_json_number(out, rows[i].x);
            out << ",\"u\":"; write_json_number(out, rows[i].u);
            out << ",\"v\":"; write_json_number(out, rows[i].v);
            out << ",\"diff\":"; write_json_number(out, rows[i].diff);
            out << '}';
        }
        out << "]}\n";
        std::cout << out.str();
    } catch (const std::exception& exc) {
        fail(exc.what());
    }
    return 0;
}
