#include <xsf/stats.h>
#include <cstddef>
struct Input {
    using value_type = double;
    const double* values;
    std::size_t n;
    std::size_t extent(int) const { return n; }
    double operator()(std::size_t i) const { return values[i]; }
};
struct Output {
    using value_type = double;
    double* values;
    std::size_t n;
    std::size_t extent(int) const { return n; }
    double& operator()(std::size_t i) const { return values[i]; }
};
extern "C" __attribute__((visibility("default")))
void run_logsf(const double* values, std::size_t n, double* result) {
    xsf::poisson_binom_logsf_all(Input{values, n}, Output{result, n + 1});
}
