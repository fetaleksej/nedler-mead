#include <iostream>
#include <vector>
#include <tuple>
#include <numeric>

#include "nedler-mead.hpp"

using namespace angstrom::math;

std::vector<std::pair<double, double>>
    generareScatteredData(double a, double b, double c)
{
    std::vector<std::pair<double, double>> data(100);
    std::accumulate(data.begin(), data.end(), 0., [a, b, c](double x, auto& item){
        double y = a * x * x + b * x + c;
        item = std::make_pair(x, y);
        return x + 0.1;
    });
    return std::move(data);
}

int main(int, char**) 
{
    auto scattered = generareScatteredData(-1, 20, 8);
    double a, b, c;
    NelderMead<3> nelder(0.,0.,0.);
    std::tie(a, b, c) = nelder.minimizing(1e-6, [scattered](double a, double b, double c) {
        return std::accumulate(scattered.begin(), scattered.end(), 0.,
                                      [a, b, c](double sum, auto& item){
            double x = item.first;
            double y = a * x * x + b * x + c;
            double e = item.second - y;
            return sum + e * e;
        });
    });

    std::cout << "Result calculation factrors "
              << "a = " << a << "; "
              << "b = " << b << "; "
              << "c = " << c << "; " << std::endl;
}
