#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

double c2 = 0.45, A = -2, B = 2, C = -2;
int callamount = 0;

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
    callamount++;
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

std::vector<double> rk2_step(double x0, const std::vector<double> &y0, double h) { //Метод из варианта
    std::vector<double> k1 = f(x0, y0);
    if (std::isnan(k1[0])) {
        return std::vector<double>(4, std::numeric_limits<double>::quiet_NaN());
    }
    std::vector<double> k2 = f(x0 + c2 * h, y0 + c2 * h * k1);
    return y0 + h * ((1 - 1 / (2 * c2)) * k1 + (1 / (2 * c2)) * k2);
}

std::vector<double> heun_step(double x0, const std::vector<double> &y0, double h) { //Метод-оппонент
    std::vector<double> k1 = f(x0, y0);
    std::vector<double> k2 = f(x0 + (1.0 / 3) * h, y0 + (1.0 / 3) * h * k1);
    std::vector<double> k3 = f(x0 + (2.0 / 3) * h, y0 + (2.0 / 3) * h * k2);
    if (std::isnan(k1[0]) or std::isnan(k2[0]) or std::isnan(k3[0])) {
        return std::vector<double>(4, std::numeric_limits<double>::quiet_NaN());
    }
    return y0 + h * (0.25 * k1 + 0.75 * k3);
}

std::vector<double> euler_step(double x0, const std::vector<double> &y0, double h) {
    return y0 + h * f(x0, y0);
}

void
auto_step_choice(std::vector<double>(*rk_step)(double, const std::vector<double> &, double), double p,
                 std::vector<double> &y0, double x0, double xk, double atol, double rtol, double h_max,
                 std::string filename1, std::string filename2, std::string filename3, std::string filename4) {
    double h, x_h2 = 0, norm, delta;
    std::vector<double> y_h2, y_h1, r_h2, r_h1;
    std::ofstream fout1(filename1), fout2(filename2), fout3(filename3), fout4(filename4);
    //Выбор начального шага
    delta = std::pow(1.0 / std::max(std::abs(x0), std::abs(xk)), p + 1) + std::pow(eunorm(f(x0, y0)), p + 1);
    h = std::pow((atol + rtol * eunorm(y0)) / delta, 1.0 / (p + 1));
    y_h1 = euler_step(x0, y0, h);
    delta = std::pow(1.0 / std::max(std::abs(x0), std::abs(xk)), p + 1) + std::pow(eunorm(f(x0 + h, y_h1)), p + 1);
    h = std::pow((atol + rtol * eunorm(y_h1)) / delta, 1.0 / (p + 1));
    fout1 << std::setprecision(16) << std::fixed << x0 << ',' << y0[0] << ',' << y0[1] << ',' << y0[2] << ',' << y0[3]
          << '\n';
    while (x0 < xk) {
        h = std::min(xk - x0, h);
        y_h1 = rk_step(x0, y0, h);
        if (std::isnan(y_h1[0])) {
            break;
        }
        y_h2 = rk_step(x0, y0, h / 2);
        x_h2 += h / 2;
        if (std::isnan(y_h2[0])) {
            break;
        }
        y_h2 = rk_step(x_h2, y_h2, h / 2);
        x_h2 += h / 2;
        double tol = (atol + rtol * eunorm(y0));
        if (std::isnan(y_h2[0])) {
            break;
        }
        r_h1 = 1.0 / (1.0 - std::pow(2.0, -p)) * (y_h2 - y_h1);
        r_h2 = 1.0 / std::pow(2.0, p) * r_h1;
        norm = eunorm(r_h2);
        if (norm > tol) {
            //fout3 << std::setprecision(16) << std::fixed << x0 << ',' << h << '\n';
            h /= 2;
        } else if ((tol / std::pow(2.0, p) < norm) and (norm <= tol)) {
            y0 = y_h2;
            x0 += h;
            //fout2 << std::setprecision(16) << std::fixed << x0 << ',' << h << '\n';
            //fout4 << std::setprecision(16) << std::fixed << x0 << ',' << eunorm(y(x0) - y0) << '\n';
            h /= 2;
            //fout1 << std::setprecision(16) << std::fixed << x0 << ',' << y0[0] << ',' << y0[1] << ',' << y0[2] << ','
                  //<< y0[3] << '\n';
        } else if ((tol / std::pow(2.0, 2 * p + 1) < norm) and (norm <= tol / std::pow(2.0, p))) {
            //fout2 << std::setprecision(16) << std::fixed << x0 << ',' << h << '\n';
            y0 = y_h1;
            x0 += h;
            //fout1 << std::setprecision(16) << std::fixed << x0 << ',' << y0[0] << ',' << y0[1] << ',' << y0[2] << ','
                  //<< y0[3] << '\n';
            //fout4 << std::setprecision(16) << std::fixed << x0 << ',' << eunorm(y(x0) - y0) << '\n';
        } else {
            //fout3 << std::setprecision(16) << std::fixed << x0 << ',' << h << '\n';
            h = std::min(2 * h, h_max);
        }
    }
}

void
optimal_step_choice(std::vector<double>(*rk_step)(double, const std::vector<double> &, double), double p,
                    std::vector<double> &y0, double x0, double xk, double tol, std::string filename) {
    std::ofstream fout(filename);
    double h, x_h1 = x0, x_h2 = x0;
    std::vector<double> y_h2 = y0, y_h1 = y0, r_h2, r_h1;
    //Выбор начального шага, для меня алгоритм выбора начального шага даёт аргумент функии, приводящий к ошибке
    /*double delta = std::pow(1.0 / std::max(std::abs(x0), std::abs(xk)), p + 1) + std::pow(eunorm(f(x0, y0)), p + 1);
    h = std::pow(tol / delta, 1.0 / (p + 1));
    y_h1 = euler_step(x0, y0, h);
    delta = std::pow(1.0 / std::max(std::abs(x0), std::abs(xk)), p + 1) + std::pow(eunorm(f(x0 + h, y_h1)), p + 1);
    h = std::pow(tol / delta, 1.0 / (p + 1));
    y_h1 = y0; //Возвращаем в исходное состояние после выбора шага */
    h = 0.015625;
    for (int i = 0; i < xk / h; i++) {
        y_h1 = rk_step(x_h1, y_h1, h);
        if (std::isnan(y_h1[0])) {
            break;
        }
        x_h1 += h;
    }
    for (int i = 0; i < xk * 2 / h; i++) {
        y_h2 = rk_step(x_h2, y_h2, h / 2);
        if (std::isnan(y_h2[0])) {
            break;
        }
        x_h2 += h / 2;
    }
    r_h2 = (1.0 / 3) * (y_h2 - y_h1);
    h = (h / 2) * sqrt(tol / eunorm(r_h2));
    while (x0 < xk) {
        h = std::min(h, xk - x0);
        y0 = rk_step(x0, y0, h);
        if (std::isnan(y0[0])) {
            break;
        }
        x0 += h;
        //yh.push_back(y0);
        fout << std::fixed << std::setprecision(16) << x0 << ',' << eunorm(y(x0) - y0) << '\n';
    }
}

int main() {
    double x0 = 0;
    std::vector<double> y0 = y(x0), y0o = y(x0), yr;
    std::vector<std::vector<double>> yh;
    std::cout.precision(16);
    std::string mode;
    std::cin >> mode;
    if (mode == "fs") { //Фиксированный шаг
        std::ofstream fout(R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output.csv)");
        for (int k = 6; k <= 12; k++) {
            double h = 1 / std::pow(2, k);
            std::cout << "Stepsize: " << std::fixed << h << std::endl;
            for (int i = 0; i <= 5.0 / h; i++) {
                y0 = rk2_step(x0, y0, h);
                y0o = heun_step(x0, y0o, h);
                if (std::isnan(y0[0]) or std::isnan(y0o[0])) {
                    break;
                }
                x0 += h;
                //yh.push_back(y0);
                yr = y(x0);
            }
            std::cout << eunorm(yr - y0) << std::endl;
            fout << std::setprecision(16) << std::fixed << h << ',' << eunorm(yr - y0) << ',' << eunorm(yr - y0o)
                 << '\n';
            y0 = std::vector<double>{1, 1, -2, 1};
            y0o = y0;
            x0 = 0;
        }
    } else if (mode == "as") { //Автоматический выбор шага с помощью оценки локальной порешности по Рунге
        auto_step_choice(rk2_step, 2, y0, 0, 5, 1e-12, 1e-6, 1,
                         R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output4.csv)",
                         R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output5.csv)",
                         R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output6.csv)",
                         R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output7.csv)");
        auto_step_choice(heun_step, 3, y0o, 0, 5, 1e-12, 1e-6, 1,
                         R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output8.csv)",
                         R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output9.csv)",
                         R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output10.csv)",
                         R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output11.csv)");
    } else if (mode == "opt") { //Вычисление оптимального шага с помощью оценки полной погрешности по Рунге
        optimal_step_choice(rk2_step, 2, y0, 0, 5, 1e-5,
                            R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output2.csv)");
        optimal_step_choice(heun_step, 3, y0o, 0, 5, 1e-5,
                            R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output3.csv)");
    } else if (mode == "as_calls") { //Самая затратная по времени операция, от 15 до 30 минут
        double rtol = 1e-4;
        std::ofstream fout1(R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output12.csv)"),
                fout2(R"(E:\All\Work\Code\SPbU\NumericMethods\DifEq\output13.csv)");
        for (int i = 0; i < 5; i++) {
            auto_step_choice(rk2_step, 2, y0, 0, 5, 1e-12, rtol, 1,
                             R"(...)", R"(...)", R"(...)", R"(...)");
            fout1<< std::setprecision(16) << std::fixed << rtol << ',' << callamount << '\n';
            callamount = 0;
            auto_step_choice(heun_step, 3, y0o, 0, 5, 1e-12, rtol, 1,
                             R"(...)", R"(...)", R"(...)", R"(...)");
            fout2<< std::setprecision(16) << std::fixed << rtol << ',' << callamount << '\n';
            callamount = 0;
            std::cout << "done with " << rtol << std::endl;
            rtol /= 10;
        }
    }
    return 0;
}
