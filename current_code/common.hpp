#pragma once

#include <cmath>
#include <vector>
#include <chrono>
#include <string>
#include <iostream>
#include <complex>

typedef std::complex<double> Complex;

// Gauss–Legendre nodes and weights on [-1, 1].
// Newton on P_n with the Numerical Recipes initial guess; fine for n ≲ 64
// (the panel code uses n_u, n_s around 8–24). Endpoint weights of a single
// n=10000 rule would need higher precision — we no longer use one.
inline void gauss_legendre(int n, std::vector<double>& x, std::vector<double>& w) {
    x.resize(n);
    w.resize(n);
    const double eps = 1.0e-15;
    const int m = (n + 1) / 2;
    for (int i = 0; i < m; ++i) {
        double z = std::cos(M_PI * (i + 0.75) / (n + 0.5));
        double z1, pp, p1;
        do {
            p1 = 1.0;
            double p2 = 0.0;
            for (int j = 0; j < n; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (std::abs(z - z1) > eps);
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
}

struct WallTimer {
    using clock = std::chrono::steady_clock;
    clock::time_point t0;
    WallTimer() : t0(clock::now()) {}
    void reset() { t0 = clock::now(); }
    double elapsed() const {
        return std::chrono::duration<double>(clock::now() - t0).count();
    }
};

inline void print_timing(const std::string& name, double seconds) {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(6);
    std::cout << "TIMING " << name << "=" << seconds << std::endl;
}

// Integer p-grid: p_i = i*dp, i = 0..np-1. Avoids the original
// `p += 0.001` loop, whose trip count depends on binary rounding of 0.001.
inline int p_count(double pmax, double dp) {
    if (dp <= 0.0) return 0;
    return static_cast<int>(std::floor(pmax / dp + 0.5)) + 1;
}
