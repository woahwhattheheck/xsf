#include <xsf/stats.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

template<class T> struct View {
    using value_type = T;
    T* data;
    std::size_t size;
    std::size_t stride = 1;
    std::size_t extent(int) const { return size; }
    T& operator()(std::size_t i) const { return data[i * stride]; }
};

template<class T>
void tables(const std::vector<T>& probabilities, int stride, std::vector<T>& pmf,
            std::vector<T>& cdf, std::vector<T>& sf, std::vector<T>& logsf) {
    const auto n = probabilities.size();
    std::vector<T> input(n * stride, T(-123));
    for (std::size_t i = 0; i < n; ++i) input[i * stride] = probabilities[i];
    View<T> p{input.data(), n, std::size_t(stride)};
    pmf.assign((n + 1) * stride, T(-123));
    cdf = sf = logsf = pmf;
    xsf::poisson_binom_pmf_all(p, View<T>{pmf.data(), n + 1, std::size_t(stride)});
    xsf::poisson_binom_cdf_all(p, View<T>{cdf.data(), n + 1, std::size_t(stride)});
    xsf::poisson_binom_sf_all(p, View<T>{sf.data(), n + 1, std::size_t(stride)});
    xsf::poisson_binom_logsf_all(p, View<T>{logsf.data(), n + 1, std::size_t(stride)});
    for (std::size_t i = 0; i < input.size(); ++i) {
        if (i % stride && input[i] != T(-123)) throw std::runtime_error("input gap changed");
    }
    for (auto* output : {&pmf, &cdf, &sf, &logsf}) {
        for (std::size_t i = 0; i < output->size(); ++i)
            if (i % stride && (*output)[i] != T(-123))
                throw std::runtime_error("output stride gap changed");
    }
}

std::vector<long double> reference_pmf(const std::vector<double>& p) {
    std::vector<long double> mass(p.size() + 1, 0);
    mass[0] = 1;
    for (std::size_t i = 0; i < p.size(); ++i) {
        long double q = 1 - static_cast<long double>(p[i]);
        for (std::size_t j = i + 1; j > 0; --j)
            mass[j] = mass[j] * q + mass[j - 1] * p[i];
        mass[0] *= q;
    }
    return mass;
}

volatile double sink = 0;
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
void calculate_logsf(View<double> p, View<double> out) {
    xsf::poisson_binom_logsf_all(p, out);
}

double median_microseconds(std::vector<double>& p, int repetitions) {
    std::vector<double> out(p.size() + 1);
    std::vector<double> measurements;
    for (int run = 0; run < 5; ++run) {
        auto start = std::chrono::steady_clock::now();
        for (int r = 0; r < repetitions; ++r) {
            calculate_logsf({p.data(), p.size()}, {out.data(), out.size()});
            sink = out[p.size() / 2];
        }
        auto end = std::chrono::steady_clock::now();
        measurements.push_back(
            std::chrono::duration<double, std::micro>(end - start).count() / repetitions
        );
    }
    std::sort(measurements.begin(), measurements.end());
    return measurements[2];
}

struct Workload {
    std::string name;
    std::vector<double> p;
    std::size_t k;
};

int main() {
    try {
        std::vector<Workload> jobs;
        std::vector<double> first;
        for (int i = 0; i < 32; ++i) first.push_back((1 + i % 9) / 128.0);
        jobs.push_back({"issue-32", first, 24});
        jobs.push_back({"issue-256", std::vector<double>(256, 1 / 128.0), 250});
        jobs.push_back({"issue-tiny-16", std::vector<double>(16, 1e-20), 1});
        std::vector<double> padded = first;
        padded.insert(padded.end(), 496, 0.0);
        padded.insert(padded.end(), 496, 1.0);
        jobs.push_back({"issue-32-with-deterministic-trials", padded, 520});
        jobs.push_back({"all-zero", std::vector<double>(1024, 0.0), 1});
        jobs.push_back({"all-one", std::vector<double>(1024, 1.0), 1024});
        jobs.push_back({"empty", {}, 1});
        std::cout << "case,n,k,sf,logsf,pmf_abs_error,cdf_abs_error,sf_abs_error,"
                     "logsf_abs_error,microseconds,stride_checked\n";
        std::cout << std::setprecision(17);
        for (auto& job : jobs) {
            std::vector<double> pmf, cdf, sf, logsf;
            tables(job.p, 3, pmf, cdf, sf, logsf);
            auto mass = reference_pmf(job.p);
            long double cumulative = 0, upper = 0;
            std::vector<long double> tail(mass.size());
            for (std::size_t j = mass.size(); j-- > 0;) {
                tail[j] = upper;
                upper += mass[j];
            }
            double pmf_error = 0, cdf_error = 0, sf_error = 0, logsf_error = 0;
            for (std::size_t j = 0; j < mass.size(); ++j) {
                if (!std::isfinite(pmf[j * 3]) || !std::isfinite(cdf[j * 3]) || !std::isfinite(sf[j * 3]))
                    throw std::runtime_error("non-finite probability table: " + job.name);
                cumulative += mass[j];
                pmf_error = std::max(pmf_error, double(std::abs(pmf[j * 3] - mass[j])));
                cdf_error = std::max(cdf_error, double(std::abs(cdf[j * 3] - cumulative)));
                sf_error = std::max(sf_error, double(std::abs(sf[j * 3] - tail[j])));
                if (tail[j] == 0) {
                    if (sf[j * 3] != 0 || logsf[j * 3] != -std::numeric_limits<double>::infinity())
                        throw std::runtime_error("wrong support boundary: " + job.name);
                } else {
                    double expected = double(std::log(tail[j]));
                    if (!std::isfinite(logsf[j * 3])) throw std::runtime_error("non-finite feasible log-tail");
                    logsf_error = std::max(logsf_error, std::abs(logsf[j * 3] - expected));
                }
            }
            if (pmf_error > 3e-14 || cdf_error > 3e-14 || sf_error > 3e-14 || logsf_error > 2e-10)
                throw std::runtime_error("reference comparison failed: " + job.name);
            const std::size_t index = job.k ? (job.k - 1) * 3 : 0;
            double elapsed = median_microseconds(job.p, job.p.size() >= 1000 ? 5 : 30);
            std::cout << job.name << ',' << job.p.size() << ',' << job.k << ','
                      << sf[index] << ',' << logsf[index] << ',' << pmf_error << ','
                      << cdf_error << ',' << sf_error << ',' << logsf_error << ','
                      << elapsed << ",3\n";
        }
        // The same original issue inputs also exercise the templated float path.
        for (const auto& job : jobs) {
            std::vector<float> p(job.p.begin(), job.p.end()), pmf, cdf, sf, logsf;
            tables(p, 2, pmf, cdf, sf, logsf);
            std::vector<double> exact_float(p.begin(), p.end());
            auto mass = reference_pmf(exact_float);
            long double upper = 0;
            for (std::size_t j = mass.size(); j-- > 0;) {
                if (upper == 0) {
                    if (sf[j * 2] != 0 || logsf[j * 2] != -std::numeric_limits<float>::infinity())
                        throw std::runtime_error("float support boundary failed");
                } else if (!std::isfinite(logsf[j * 2]) ||
                           std::abs(double(logsf[j * 2]) - double(std::log(upper))) > 0.03) {
                    throw std::runtime_error("float log-tail reference failed");
                }
                if (!std::isfinite(pmf[j * 2]) || !std::isfinite(cdf[j * 2]) || !std::isfinite(sf[j * 2]))
                    throw std::runtime_error("non-finite float probability table");
                upper += mass[j];
            }
        }
        std::cerr << "Original cases, all four table functions, float/double and strided views completed.\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
