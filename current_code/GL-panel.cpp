// Composite Gauss–Legendre SFA2 wavepacket integrator.
//
// Same physical integrand as original_code/GL-test.cpp, but:
//   - A, S1=∫A, S2=∫A², and k are evaluated from closed forms (no m_sub
//     trapezoid, no O(n³) k table, no bicubic interpolant of k);
//   - the (u,v) triangle is tiled with small GL panels instead of one
//     global n×n product rule;
//   - the first v-panel uses v = s² so the model-1 1/sqrt(v) end is cancelled;
//   - p-independent work is done once per node; the inner loop is over p.
//
// Coordinates: t = u+v, t' = u-v, v >= 0. The original extra factor of 2
// in front of the integrand is kept (w2 = 2 * GL_weight).
//
// Defaults n_u=8, n_s=12, n_s_panels=40, osc_per_panel=2 were checked against
// a 24/20/80/osc=1 reference: abs/peak ~ 1e-8 on p < 2.3.

#include "common.hpp"

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <omp.h>

using namespace std;

double F0, omega, pi;
int nc, model;
double alpha, beta_par, epsilon;
double period, tlength, Ip, Ip_ion;
double pmax, t0, tmax;
double alpha2, beta2;

int n_u = 8;           // GL order per u-panel
int n_s = 12;          // GL order per v-panel (or s-panel on the first slice)
int n_s_panels = 40;   // equal partitions of [0, vmax]
double osc_per_panel = 2.0;  // target Volkov oscillations per u-panel at pmax
int use_s2 = 1;        // 1: map only the first v-panel with v=s^2

// Clamp times used for A, F, S1, S2, k. Ip*t' and 1/sqrt(|t-t'|) stay unclamped,
// matching original_code/GL-test.cpp.
inline double clamp_time(double t) {
    if (t < 0.0) return 0.0;
    if (t > tmax) return tmax;
    return t;
}

// Vector potential of a sin^2 envelope × sin(omega t). Zero after the pulse.
inline double A_func(double t_val) {
    if (t_val > tlength || t_val < 0.0) return 0.0;
    const double senv = std::sin((pi * t_val) / tlength);
    return (F0 * senv * senv * std::sin(t_val * omega)) / omega;
}

// F = -dA/dt, analytic (original used this too; older SFA codes used a finite difference).
inline double Field_func(double t_val) {
    if (t_val > tlength || t_val < 0.0) return 0.0;
    const double senv = std::sin((pi * t_val) / tlength);
    const double cenv = std::cos((pi * t_val) / tlength);
    const double sosc = std::sin(t_val * omega);
    const double cosc = std::cos(t_val * omega);
    const double d_env_sq_dt = (2.0 * pi / tlength) * senv * cenv;
    return -(F0 / omega) * (d_env_sq_dt * sosc + senv * senv * omega * cosc);
}

// Primitives used by S1, S2. A is a sum of three sines after expanding sin^2.
inline double integ_sin(double w, double t) { return (1.0 - std::cos(w * t)) / w; }
inline double integ_sin2(double w, double t) { return 0.5 * t - std::sin(2.0 * w * t) / (4.0 * w); }
inline double integ_sinsin(double a, double b, double t) {
    if (std::abs(a - b) < 1.0e-14) return integ_sin2(a, t);
    return 0.5 * (std::sin((a - b) * t) / (a - b) - std::sin((a + b) * t) / (a + b));
}

// S1(t) = ∫_0^{min(t,T)} A(τ) dτ. After the pulse A=0 so S1 freezes at S1(T).
// For an integer number of cycles S1(T)=0 exactly (net vector potential vanishes).
inline double S1_analytic(double t) {
    if (t <= 0.0) return 0.0;
    const double tc = (t > tlength) ? tlength : t;
    const double wp = omega + 2.0 * pi / tlength;
    const double wm = omega - 2.0 * pi / tlength;
    const double C = F0 / (2.0 * omega);
    return C * (integ_sin(omega, tc) - 0.5 * integ_sin(wp, tc) - 0.5 * integ_sin(wm, tc));
}

// S2(t) = ∫_0^{min(t,T)} A(τ)² dτ. Expand A as three sines, integrate sin² and sin·sin.
inline double S2_analytic(double t) {
    if (t <= 0.0) return 0.0;
    const double tc = (t > tlength) ? tlength : t;
    const double wp = omega + 2.0 * pi / tlength;
    const double wm = omega - 2.0 * pi / tlength;
    const double C = F0 / (2.0 * omega);
    const double val = integ_sin2(omega, tc)
        + 0.25 * integ_sin2(wp, tc)
        + 0.25 * integ_sin2(wm, tc)
        - integ_sinsin(omega, wp, tc)
        - integ_sinsin(omega, wm, tc)
        + 0.5 * integ_sinsin(wp, wm, tc);
    return (C * C) * val;
}

// p-independent payload of one quadrature node. The p-loop only needs these.
struct NodeCache {
    double t;     // clamped t, used in 0.5 p² t
    double S1t;
    double k;     // Morishita momentum: -(S1(t)-S1(t')) / (|t-t'|+eps)
    double w2;    // 2 * (GL weight including Jacobians); the 2 is from the original integrand
    double Gpre;  // G without the (k-p)² factor
    double phi0;  // phase without 0.5 p² t + p S1(t)
};

inline void fill_node(double u, double v, double weight, NodeCache& nd) {
    const double t = u + v;
    const double tp = u - v;
    const double t_c = clamp_time(t);
    const double tp_c = clamp_time(tp);

    nd.t = t_c;
    nd.S1t = S1_analytic(t_c);
    const double S1tp = S1_analytic(tp_c);
    const double S2tp = S2_analytic(tp_c);
    const double A_tp = A_func(tp_c);
    const double F_tp = Field_func(tp_c);

    // Same k the original n×n table stored: average of -A between t' and t.
    nd.k = -(nd.S1t - S1tp) / (std::abs(t_c - tp_c) + epsilon);

    // Action(t,p)+Action(t',k)-Action(t,k)+Ip*t' with S2(t) cancelled.
    // Ip*tp uses the raw (unclamped) t', as in GL-test.cpp.
    nd.phi0 = nd.k * (S1tp - nd.S1t) + 0.5 * S2tp
            + 0.5 * nd.k * nd.k * (tp_c - t_c) + Ip * tp;
    nd.w2 = 2.0 * weight;

    const double kA = nd.k + A_tp;
    const double den_beta = beta2 + kA * kA;
    if (model == 1) {
        // sqrt uses raw |t-tp|=2v, not the clamped times.
        nd.Gpre = F_tp * kA / std::sqrt(std::abs(t - tp) + epsilon) / (den_beta * den_beta);
    } else {
        const Complex Z_aux(epsilon, -(t - tp));
        const double zabs = std::abs(Z_aux);
        const double znorm = std::norm(Z_aux);
        const double ratio = (2.0 * pi * zabs) / znorm;
        const double pref = ratio * std::sqrt(ratio); // ( )^{1.5}
        nd.Gpre = pref * F_tp * kA / (den_beta * den_beta * den_beta);
        nd.phi0 += 1.5 * std::arg(Z_aux);
    }
}

struct SNode {
    double v;
    double wdv; // weight already including ds or dv and, for s², the 2s Jacobian
};

// Build 1-D quadrature nodes in v.
// Mapping the whole [0,vmax] with v=s² and grading s as (j/M)^2 starves the
// bulk of v (that experiment was ~10–100% wrong). s² is used only on the
// first panel, where 1/sqrt(v) actually sits.
void build_s_nodes(std::vector<SNode>& snodes, const std::vector<double>& xs,
                   const std::vector<double>& ws) {
    const double vmax = 0.5 * (tmax - t0);
    snodes.clear();
    snodes.reserve(n_s_panels * n_s);
    for (int p = 0; p < n_s_panels; ++p) {
        const double va = vmax * (double)p / n_s_panels;
        const double vb = vmax * (double)(p + 1) / n_s_panels;
        if (use_s2 && p == 0) {
            const double sa = std::sqrt(std::max(va, 0.0));
            const double sb = std::sqrt(std::max(vb, 0.0));
            const double smid = 0.5 * (sa + sb);
            const double shalf = 0.5 * (sb - sa);
            for (int i = 0; i < n_s; ++i) {
                const double s = smid + shalf * xs[i];
                if (s < 0.0) continue; // ξ=-1 on [0,sb] gives s=0, not negative
                SNode nd;
                nd.v = s * s;
                nd.wdv = ws[i] * shalf * (2.0 * s); // dv = 2s ds
                snodes.push_back(nd);
            }
        } else {
            const double vmid = 0.5 * (va + vb);
            const double vhalf = 0.5 * (vb - va);
            for (int i = 0; i < n_s; ++i) {
                SNode nd;
                nd.v = vmid + vhalf * xs[i];
                nd.wdv = ws[i] * vhalf;
                if (nd.v < 0.0) nd.v = 0.0;
                snodes.push_back(nd);
            }
        }
    }
}

void integrate_all_p(const std::vector<double>& p_vals, std::vector<Complex>& results) {
    const int np = static_cast<int>(p_vals.size());
    std::vector<double> xu, wu, xs, ws;
    gauss_legendre(n_u, xu, wu);
    gauss_legendre(n_s, xs, ws);

    std::vector<SNode> snodes;
    build_s_nodes(snodes, xs, ws);

    // Panel length in u is set by the fastest phase on the scan, Ω = pmax²/2 + Ip.
    // Small-p points are therefore over-resolved; one node set is reused for all p.
    const double Omega = 0.5 * pmax * pmax + Ip;
    const double lambda = (Omega > 1e-12) ? (2.0 * pi / Omega) : (tmax - t0);
    const double u_len = std::max(lambda * osc_per_panel, 1e-6);

    long long n_s_nodes = (long long)snodes.size();
    long long n_u_est = 0;
    if (!snodes.empty()) {
        const double span = (tmax - t0);
        n_u_est = (long long)std::max(1, (int)std::ceil(span / u_len)) * n_u;
    }
    std::cout << "QUAD n_u=" << n_u << " n_s=" << n_s
              << " n_s_panels=" << n_s_panels
              << " s_nodes=" << n_s_nodes
              << " u_panel_len=" << u_len
              << " osc_per_panel=" << osc_per_panel
              << " Omega_max=" << Omega << std::endl;
    std::cout << "QUAD approx_u_points_at_v0=" << n_u_est
              << " approx_nodes=" << (n_s_nodes * n_u_est) << std::endl;
    {
        const double A0 = (omega > 0.0) ? (F0 / omega) : 0.0;
        const double Omega_kin = 0.5 * (pmax + A0) * (pmax + A0) + Ip;
        std::cout << "QUAD A0=" << A0 << " Omega_p=" << Omega
                  << " Omega_kin=" << Omega_kin
                  << " Omega_kin/Omega_p=" << (Omega > 0.0 ? Omega_kin / Omega : 0.0)
                  << std::endl;
    }

    const int ns = (int)snodes.size();
    #pragma omp parallel
    {
        // Private spectrum per thread; adding an np-array with OpenMP reduction
        // on complex<> is awkward, so we combine in a short critical section.
        std::vector<Complex> local(np, Complex(0.0, 0.0));
        #pragma omp for schedule(dynamic, 1)
        for (int is = 0; is < ns; ++is) {
            const double v = snodes[is].v;
            const double wdv = snodes[is].wdv;
            const double u_low = t0 + v;
            const double u_high = tmax - v;
            if (u_high <= u_low) continue; // v near vmax: the triangle has no width

            const int npan = std::max(1, (int)std::ceil((u_high - u_low) / u_len));
            const double h = (u_high - u_low) / npan; // tile exactly, no leftover stub
            for (int ipan = 0; ipan < npan; ++ipan) {
                const double ua = u_low + ipan * h;
                const double ub = ua + h;
                const double umid = 0.5 * (ua + ub);
                const double uhalf = 0.5 * (ub - ua);
                for (int iu = 0; iu < n_u; ++iu) {
                    const double u = umid + uhalf * xu[iu];
                    const double weight = wu[iu] * uhalf * wdv;
                    NodeCache nd;
                    fill_node(u, v, weight, nd);
                    for (int ip = 0; ip < np; ++ip) {
                        const double p = p_vals[ip];
                        const double kp = nd.k - p;
                        const double G = nd.Gpre / (alpha2 + kp * kp);
                        const double phi = nd.phi0 + (0.5 * nd.t) * p * p + nd.S1t * p;
                        const double c = std::cos(phi);
                        const double s = std::sin(phi);
                        local[ip] += Complex(nd.w2 * G * c, nd.w2 * G * s);
                    }
                }
            }
        }
        #pragma omp critical
        {
            for (int ip = 0; ip < np; ++ip) results[ip] += local[ip];
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc < 8) {
        cerr << "Usage: ./GL-panel Ip_ion nc omega F0 model pmax dp [outfile]\n"
             << "       [n_u n_s n_s_panels osc_per_panel use_s2]\n";
        return 1;
    }
    Ip_ion = std::atof(argv[1]); // read for CLI compatibility; the phase uses Ip = 0.5
    nc = std::atoi(argv[2]);
    omega = std::atof(argv[3]);
    F0 = std::atof(argv[4]);
    model = std::atoi(argv[5]);
    pmax = std::atof(argv[6]);
    const double dp = std::atof(argv[7]);
    std::string outfile = (argc >= 9) ? argv[8] : "";
    if (argc >= 10) n_u = std::atoi(argv[9]);
    if (argc >= 11) n_s = std::atoi(argv[10]);
    const bool nsp_from_cli = (argc >= 12);
    if (nsp_from_cli) n_s_panels = std::atoi(argv[11]);
    if (argc >= 13) osc_per_panel = std::atof(argv[12]);
    if (argc >= 14) use_s2 = std::atoi(argv[13]);

    Ip = 0.5; // bound-electron Ip in the Volkov phase; not Ip_ion
    pi = M_PI;
    if (model == 1) {
        alpha = 1.55; beta_par = 0.75; epsilon = 1e-8;
    } else if (model == 2) {
        alpha = 1.0; beta_par = 1.0; epsilon = 1e-4;
    } else {
        cerr << "Unknown model\n";
        return 1;
    }
    alpha2 = alpha * alpha;
    beta2 = beta_par * beta_par;
    period = 2.0 * pi / omega;
    tlength = nc * period;
    tmax = tlength * 1.5; // tail after the pulse so recollision remains possible
    t0 = 0.0;

    // v-structure of k ~ A ~ F0/ω. Calibrated at F0=ω=0.06 (A0=1 → 40 panels).
    // u-panels already follow pmax via Ω = pmax²/2+Ip; they do not need extra F.
    if (!nsp_from_cli) {
        const double A0 = (omega > 0.0) ? (F0 / omega) : 1.0;
        n_s_panels = std::max(24, (int)std::lround(40.0 * A0));
    }

    const int np = p_count(pmax, dp);
    std::cout << "PROG panel n_u=" << n_u << " n_s=" << n_s
              << " n_s_panels=" << n_s_panels << " np=" << np
              << " model=" << model << " threads=" << omp_get_max_threads() << std::endl;

    WallTimer total;
    vector<double> p_vals(np);
    for (int ip = 0; ip < np; ++ip) p_vals[ip] = ip * dp;
    vector<Complex> results(np, Complex(0.0, 0.0));

    WallTimer t;
    integrate_all_p(p_vals, results);
    print_timing("integrate_sec", t.elapsed());

    if (outfile.empty()) {
        ostringstream os;
        os << "panel_np=" << np << ".dat";
        outfile = os.str();
    }
    {
        ofstream out(outfile, ios::trunc);
        out.setf(ios::scientific);
        out.precision(16);
        // First column is p*0.55, same plotting convention as GL-test.cpp.
        for (int ip = 0; ip < np; ++ip)
            out << p_vals[ip] * 0.55 << "    " << std::norm(results[ip]) << "\n";
    }
    print_timing("total_sec", total.elapsed());
    std::cout << "OUTFILE " << outfile << std::endl;
    return 0;
}
