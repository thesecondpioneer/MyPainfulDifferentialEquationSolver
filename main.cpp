#include <iostream>
#include <cmath>
#include <vector>

double c2 = 0.45, A = -2, B = 2, C = -2;

std::vector<double> operator+(std::vector<double> a, std::vector<double> b) {
    for (int i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
    return a;
}

std::vector<double> operator-(std::vector<double> a, std::vector<double> b) {
    for (int i = 0; i < a.size(); i++) {
        a[i] -= b[i];
    }
    return a;
}

std::vector<double> operator*(double a, std::vector<double> b) {
    for (int i = 0; i < b.size(); i++) {
        b[i] *= a;
    }
    return b;
}

double eunorm(std::vector<double> a) {
    double sum(0);
    for (int i = 0; i < a.size(); i++) {
        sum += a[i] * a[i];
    }
    return std::sqrt(sum);
}

std::vector<double> f(double x, std::vector<double> y) {
    if ((y[1] < 0) or (y[0] < 0)) {
        std::cout << "Illegal argument, can't use this stepsize" << std::endl;
        return std::vector<double>(4, std::numeric_limits<double>::quiet_NaN());
    }
    return std::vector<double>{2.0 * x * std::pow(y[1], 1.0 / B) * y[3],
                               2.0 * B * x * std::exp(B / C * (y[2] - A)) * y[3],
                               2.0 * C * x * y[3],
                               -2.0 * x * std::log(y[0])
    };
}

std::vector<double> y(double x) {
    double x2 = std::pow(x, 2), sx2 = std::sin(x2), bsx2 = B * sx2, csx2 = C * sx2;
    return std::vector<double>{std::exp(sx2),
                               std::exp(bsx2),
                               csx2 + A,
                               std::cos(x2),
    };
}

std::vector<double> rk2_step(double x0, const std::vector<double> &y0, double h) {
    std::vector<double> k1 = f(x0, y0);
    if (std::isnan(k1[0])) {
        return std::vector<double>(4, std::numeric_limits<double>::quiet_NaN());
    }
    std::vector<double> k2 = f(x0 + c2 * h, y0 + c2 * h * k1);
    return y0 + h * ((1 - 1 / (2 * c2)) * k1 + (1 / (2 * c2)) * k2);
}

std::vector<double> euler_step(double x0, const std::vector<double> &y0, double h) {
    return y0 + h * f(x0, y0);
}

int main() {
    double x0 = 0, xk = 5.0;
    std::vector<double> y0 = y(x0), yr;
    std::vector<std::vector<double>> yh;
    std::cout.precision(16);
    std::string mode;
    std::cin >> mode;
    if (mode == "fs") {
        for (int k = 1; k <= 6; k++) {
            double h = 1 / std::pow(2, k);
            std::cout << "Stepsize: " << std::fixed << h << std::endl;
            for (int i = 0; i <= 5.0 / h; i++) {
                y0 = rk2_step(x0, y0, h);
                if (std::isnan(y0[0])) {
                    break;
                }
                x0 += h;
                yh.push_back(y0);
                yr = y(x0);
                std::cout << eunorm(yr - y0) << std::endl;
            }
            y0 = std::vector<double>{1, 1, -2, 1};
            x0 = 0;
        }
    } else if (mode == "as") {
        const double h_max = 1.0;
        double h, x_h2 = 0, norm, delta, atol = 1e-12, rtol = 1e-8;
        std::vector<double> y_h2, y_h1, r_h2, r_h1;
        delta = std::pow(1.0 / std::max(std::abs(x0), std::abs(xk)), 3) + std::pow(eunorm(f(x0, y0)), 3);
        h = std::pow((atol + rtol * eunorm(y0)) / delta, 1.0 / 3.0);
        y_h1 = euler_step(x0, y0, h);
        delta = std::pow(1.0 / std::max(std::abs(x0), std::abs(xk)), 3) + std::pow(eunorm(f(x0 + h, y_h1)), 3);
        h = std::pow((atol + rtol * eunorm(y_h1)) / delta, 1.0 / 3.0);
        while (x0 < 5.0) {
            h = std::min(5.0 - x0, h);
            y_h1 = rk2_step(x0, y0, h);
            if (std::isnan(y_h1[0])) {
                break;
            }
            y_h2 = rk2_step(x0, y0, h / 2);
            x_h2 += h / 2;
            if (std::isnan(y_h2[0])) {
                break;
            }
            y_h2 = rk2_step(x_h2, y_h2, h / 2);
            x_h2 += h / 2;
            double tol = (atol + rtol * eunorm(y0));
            if (std::isnan(y_h2[0])) {
                break;
            }
            r_h1 = (4.0 / 3.0) * (y_h2 - y_h1);
            r_h2 = (1.0 / 4.0) * (r_h1);
            norm = eunorm(r_h2);
            if (norm > tol) {
                h /= 2;
                continue;
            } else if ((tol / 4 < norm) and (norm <= tol)) {
                y0 = y_h2;
                x0 += h;
                h /= 2;
                continue;
            } else if ((tol / 32 < norm) and (norm <= tol / 4)) {
                y0 = y_h1;
                x0 += h;
                continue;
            } else {
                h = std::min(2 * h, h_max);
                continue;
            }
        }
        std::cout << eunorm(y(5.0) - y0) << std::endl;
    } else if (mode == "opt") {
        double h = 0.015625, x_h1 = x0, x_h2 = x0, h_prev = INT64_MAX, tol = 0.00001;
        std::vector<double> y_h2 = y0, y_h1 = y0, r_h2, r_h1;
        //std::cin >> tol;
        for (int i = 0; i < 5.0 / h; i++) {
            y_h1 = rk2_step(x_h1, y_h1, h);
            if (std::isnan(y_h1[0])) {
                break;
            }
            x_h1 += h;
        }
        for (int i = 0; i < 10.0 / h; i++) {
            y_h2 = rk2_step(x_h2, y_h2, h / 2);
            if (std::isnan(y_h2[0])) {
                break;
            }
            x_h2 += h / 2;
        }
        r_h2 = (1.0 / 3) * (y_h2 - y_h1);
        h = (h / 2) * sqrt(tol / eunorm(r_h2));
        for (int i = 0; i <= 5.0 / h; i++) {
            y0 = rk2_step(x0, y0, h);
            if (std::isnan(y0[0])) {
                break;
            }
            x0 += h;
            yh.push_back(y0);
            yr = y(x0);
        }
        std::cout << eunorm(yr - y0) << std::endl;
    }
    return 0;
}
