# Nedler-Mead
It is project have been created for demonstration processing Nelder-Mead algorithm
Following example shows working nedler-mead algorithm:
```c++
int main(int, char**) 
{
    // In result we get vector of pairs where first is 'x' and second is 'y'
    auto scattered = generareScatteredData(-1, 20, 8);
    double a, b, c;
    NelderMead<3> nelder(0.,0.,0.); // Create worker object
    // At bellow line we call function 'minimizing' with lambda argument
    // where we search minimal value of lambda
    std::tie(a, b, c) = nelder.minimizing(1e-6, [scattered](double a, double b, double c) {
        // We calculate deviation between current function value and expected value
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

```

