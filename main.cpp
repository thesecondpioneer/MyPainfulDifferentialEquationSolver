#include <iostream>
#include <cmath>
#include <vector>

double c2 = 0.45, A = -2, B = 2, C = -2;

std::vector<double> operator+(std::vector<double> a, std::vector<double> b){
    for(int i = 0; i < a.size(); i++){
        a[i] += b[i];
    }
    return a;
}

std::vector<double> operator-(std::vector<double> a, std::vector<double> b){
    for(int i = 0; i < a.size(); i++){
        a[i] -= b[i];
    }
    return a;
}

std::vector<double> operator*(double a, std::vector<double> b){
    for(double i:b){
        i *= a;
    }
    return b;
}

double eunorm(std::vector<double> a){
    double sum(0);
    for (int i = 0; i < a.size(); i++){
        sum += a[i] * a[i];
    }
    return std::sqrt(sum);
}

std::vector<double> f(double x, std::vector<double> y) {
    return std::vector<double>{2.0 * x * std::pow(std::abs(y[1]), 1.0 / B) * y[3],
                               2.0 * B * x * std::exp(B / C * (y[2] - A)) * y[3],
                               2.0 * C * x * y[3],
                               -2.0 * x * std::log(std::abs(y[0])),
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

std::vector<double> rk2_step(double x0, std::vector<double> y0, double h){
    std::vector<double> k1 = f(x0, y0);
    return y0 + h * ((1 - 1 / (2 * c2)) * k1 + (1 / (2 * c2)) * f(x0 + c2 * h, y0 + c2 * h * k1));
}

std::vector<double> euler_step(double x0, std::vector<double> y0, double h){
    return y0 + h * f(x0, y0);
}

int main() {
    double x0 = 0;
    std::vector<double> y0 = y(x0), yr;
    std::vector<std::vector<double>> yh;
    std::cout.precision(16);
    for (int k = 11; k <= 11; k++){
        double h = 1/std::pow(2, k);
        std::cout << "Stepsize: " << std::fixed << h << std::endl;
        for(int i = 1; i <= 5.0 / h; i++) {
            y0 = rk2_step(x0, y0, h);
            x0 += h;
            yh.push_back(y0);
            yr = y(x0);
            std::cout << eunorm(yr - y0) << std::endl;
        }
        y0 = std::vector<double>{1, 1, -2, 1};
        x0 = 0;
    }
    return 0;
}
