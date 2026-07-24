/*************************************************************
* SFA2 solver for obtaining returning electronic wavepackets *
* author: Dmitry Efimov 2024,2025                            *
**************************************************************/
//wczytujemy biblioteki
#include <iostream>
#include <complex>
#include <functional>
#include <math.h>
#include <vector>
#include <ctime>
#include <complex>
#include <cmath>
#include <fstream>  // Ensure this is included
#include <omp.h>
#include <mutex>  // Ensure this is included
#include <fstream>

using namespace std;

const complex<double> ImagUnity(0, 1);
const complex<double> Iu(0, 1);
const    double Pi = M_PI;
const double TOL   = 1e-6; // threshold

typedef std::complex<double> Complex;

// Declare the reduction operation for std::complex<double>
#pragma omp declare reduction(+: Complex: omp_out += omp_in) initializer(omp_priv = Complex(0.0, 0.0))
//sumowanie liczb zespolonych
std::mutex print_mutex;
std::complex<double> im_un(0.0, 1.0); // imaginary unity
// deklarujemy zmienne
    double F0,xmin,xmax,ymin,ymax;
    double omega;
    double pi;
    int nc,n,N,model;
    double alpha;
    double beta_par;
    double epsilon;
    double p_current;
    double period;
    double tlength;
    double Ip;
    double Ip_ion;
    double PI;
    double F0_nominal;
    double pmax;
    int zero_above_diagonal;
    double dt;
    int test_of_phys;
    double t0,tmax;

std::vector<double> A_array;
std::vector<double> F_array;
std::vector<double> S1_array;
std::vector<double> S2_array;
std::vector<double> k_array;
std::vector<double> check_array;


void exportVectorToFile(const std::vector<double>& vec, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    for (const auto& val : vec) {
        out << val << "\n";
    }
}

void export2DArray(const std::vector<double>& data, int rows, int cols, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int idx = cols * i + j;
            out << j << " " << i << " " << data[idx] << "\n";
        }
        out << "\n"; // blank line for gnuplot's `splot` with grid mode
    }
}

void writeNumbersToFile(const std::string &filename, double num1, double num2) {
    // Create an output file stream (ofstream) object in append mode
    std::ofstream outFile;

    // Open the file in append mode
    outFile.open(filename, std::ios::app);

    // Check if the file is open
    if (!outFile.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return;
    }

    // Write the numbers to the file
    outFile << num1 << "    " << num2 << "    " << std::endl;

    // Close the file
    outFile.close();

    std::cout << "Numbers written to file successfully for p = " << num1 << " that is " << num2 << std::endl;
}

void writeNumbersToFileErr(const std::string &filename, double num1) {
    // Create an output file stream (ofstream) object in append mode
    std::ofstream outFile;

    // Open the file in append mode
    outFile.open(filename, std::ios::app);

    // Check if the file is open
    if (!outFile.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return;
    }

    // Write the numbers to the file
    outFile << num1  << std::endl;

    // Close the file
    outFile.close();

    //std::cout << "Numbers written to file successfully!" << std::endl;
}



//****************************************************************************************************************************************
//****** PREPARATION OF THE FUNCTION TO BE INTEGRATED ************************************************************************************
//****************************************************************************************************************************************


int isPositive(double number) {
    return number > 0.0;
}
//definiujemy funkcje signum
inline double Sign(double x) {
    if (x > 0) return 1.;
    if (x < 0) return -1.;
    return 0.; // Return 0 if x is zero
}

inline int index2D(int i, int j){
    return n*i + j;
}
//potencjal wektorowy
// vector potential
double A(double t, double F0, double omega, double tlength) {
    if (t > tlength) {
        return 0.0;
    }
    else {
        double T=tlength;
        return (F0 * std::pow(std::sin((pi * t) / T), 2.0) * std::sin(t * omega)) / omega;
    }
}

// field value
//pochodna po potencjale wektorowym- indukowane pole
double Field(double t, double F0, double omega, double tlength) {
    if (t > tlength) {
        return 0.0;
    }
    else { //definiujemy mały przyrost różniczke dtt
        double dtt = 0.0000001;
        return -(A(t+dtt, F0, omega, tlength) - A(t, F0, omega, tlength))/dtt;
    }
}
//obliczanie potencjału wektorowego?
// computing array of vector potential values for the nodes of computation grid
void compute_A_array(){
    for (int i = 0; i < n; i += 1){
        A_array[i]= A(t0 + i*dt, F0, omega, tlength);
    }
    cout << "A done" << endl;
}

// computing array of electric field values for the nodes of computation grid
void compute_F_array(){
    for (int i = 0; i < n; i += 1){
        F_array[i]= Field(i*dt, F0, omega, tlength);
    }
    cout << "F done" << endl;
}

// S1_array is an integral from 0 to i*dt of A(t) - vector potential
// k_array stands for values of saddle points  - see Efimov,2025 or Morishita et al 2007
void compute_S1_and_k_array(){
    int m = 10000; // a factor to wich we divide an elementary square on t,tp plane
    double dt_loc = dt/(m-1);
    double integ;
    for (int k = 0; k < n; k += 1){
        S1_array[k] = 0.0;
        for (int k2 = 0; k2 < n; k2 += 1){
            k_array[index2D(k2,k)] = 0.0;
        }
    }
    // for each elementary time interval we compute integral of vector potential over it
    for (int i = 1; i < n; i += 1){
        integ = -(A((i-1)*dt, F0, omega, tlength) + A(i*dt, F0, omega, tlength))/2.0;
        for (int j = 0; j < m; j += 1){
            integ += A(((i-1)*dt + j*dt_loc), F0, omega, tlength);
        }
        integ *= dt_loc;
        for (int kx = 0; kx < n; kx += 1){
            // if we find that this elementary interval is within S1 integration interval (0,t), we add the computed elementary integral value to S1 integral
            if (i <= kx){
                S1_array[kx] += integ;
            }
            for (int ky = 0; ky < n; ky += 1){
                // if we find that this elementary interval is within k_array integration interval, that is from tp to t (tp <-> y , t <-> x)
                // we add the computed elementary integral value to k_array integral
                if ((i-1) >= ky and i <= kx and ky < kx){
                    k_array[index2D(kx,ky)] += integ;
                }
            }
        }
        /*
        // screening y or tp
        for (int k = i; k < n; k += 1){
            S1_array[k] += integ;
            //screening x or t
            //for (int k2 = 0; k2 < i; k2 += 1){
            for (int k2 = 1; k2 <= i; k2 += 1){
                k_array[index2D(k2,k-1)] += integ;
                //k_array[index2D(k-1,k2)] += integ;
            }
        }
        */
    }

    for (int k = 0; k < n; k += 1){
        for (int k2 = 0; k2 < n; k2 += 1){
            k_array[index2D(k2,k)] *= -1.0/ (abs(k2*dt - k*dt) + epsilon); //(abs(t - tp) + epsilon);
        }
    }

    cout << "S1 done" << endl;
}
//obliczamy działanie?
// S2_array is an integral from 0 to i*dt of A^2(t) - square of vector potential
void compute_S2_array(){
    int m = 10000;
    double dt_loc = xmax/(n-1)/(m-1);
    double integ;
    for (int k = 0; k < n; k += 1){
        S2_array[k] = 0.0;
    }
    for (int i = 1; i < n; i += 1){
        integ = -(A((i-1)*dt, F0, omega, tlength) * A((i-1)*dt, F0, omega, tlength) + A(i*dt, F0, omega, tlength) * A(i*dt, F0, omega, tlength))/2.0;
        for (int j = 0; j < m; j += 1){
            integ += A(((i-1)*dt + j*dt_loc), F0, omega, tlength) * A(((i-1)*dt + j*dt_loc), F0, omega, tlength);
        }
        integ *= dt_loc;
        for (int k = i; k < n; k += 1){
            S2_array[k] += integ;
        }
    }
    cout << "S2 done" << endl;
}

// Compute action as a function of time (descrete parameter i: t = i*dt) and momentum p.
// ! The integration is from 0 to a currect time passed as a parameter
//definiujemy działanie?
double Action(int i, double p) {
    double int1 = 0.0;
    double int2 = 0.0;
    double int3 = 0.0;

    int1 = p*p*i*dt/2.0;
    //int1 = i*dt/2.0;

    int2 = p * S1_array[i];

    int3 = 1.0/2.0 * S2_array[i];

    // === Final result ===
    return int1 + int2 + int3;
    //return int3;// + int2 + int3;
}
// elektron znajduje się w określonych fazach?
// Compute the phase function for given parameters
double Phase_saddle(double p, double Ip, double tp, double k, int i, int j) {
    return Action(i, p) + Action(j, k) - Action(i, k) + Ip * tp;
}

void checking_action(){
    double k;
    for (int i = 0; i < n; i += 1){
        for (int j = 0; j < n; j += 1){
            k = k_array[index2D(i,j)];
            check_array[index2D(i,j)] = Phase_saddle(0.0, Ip, j*dt, k, i, j);
        }
    }

    cout << "checking x done" << endl;
}

// computation of amplitude and phase of a SFA2 expression
std::vector<double> generatePoints_phys(int x, int y) {
    double p=p_current;
    double k;
    double rho = 0.5;
    double a_prm = 1.0;
    double t = t0 + x * dt;
    double tp = t0 + y * dt;

    int ii = x;
    int jj = y;
    double G_here;
    double Phi_here;
    double Z_charge = beta_par;
    complex<double>  Z_aux=0.0;
    //k = k_array[index2D(jj,ii)] ;
    k = k_array[index2D(ii,jj)] ;
    if (model == 1){
        G_here = isPositive(t - tp) * 1./sqrt(abs(t - tp) + epsilon) * F_array[jj] * (k + A_array[jj]) / (pow(alpha, 2) + pow(k - p, 2)) / pow(pow(beta_par, 2) + pow(k + A_array[jj],2), 2);
        Phi_here = Phase_saddle(p, Ip, tp, k, ii, jj);
    }
    //!
    else if (model==2){
        Z_aux = epsilon - im_un*(t-tp);
        G_here = isPositive(t - tp) * pow((2.0*PI*abs(Z_aux))/(abs(Z_aux)*abs(Z_aux)),3.0/2.0) * F_array[jj] * (k + A_array[jj]) / (pow(alpha, 2) + pow(k - p, 2)) / pow(pow(Z_charge, 2) + pow(k + A_array[jj],2), 3);
        Phi_here = Phase_saddle(p, Ip, tp, k, ii, jj) + 3.0/2.0*arg(Z_aux);
    }
    else {
        std::cout << "FATAL ERROR: MODEL UNIDENTIFIED (FUNCTION GENERATEPOINTS)" << std::endl;
        exit(1);
    }
    std::vector<double> Vector(2);
    Vector[0] = G_here;
    Vector[1] = Phi_here;
    //std::cout << "k: " << k << ", alpha: " << alpha << ", p: " << p << " G_here  " << G_here << "  Phi_here  " << Phi_here << std::endl;
    return Vector;
}

// computation of amplitude and phase of an arbitrary test function
std::vector<double> generatePoints_test(int x, int y) {

    double t = t0 + x * dt;
    double tp = t0 + y * dt;
    double G_here;
    double Phi_here;

       // G_here = 1.0;//cos(t*tp*tp/1.);
        //Phi_here = t*tp*tp*tp*tp;
        //G_here = -1.0*(t*t+tp*tp)/100.;//cos(t*tp*tp/1.);
        //Phi_here = t*t * tp*tp;
        //G_here = 1.0;
        //Phi_here = 100*(t*t*t*2*tp*tp);
        //G_here   = 1.2 + 0.3*t - 0.7*tp;     // exactly linear
        //G_here   = -(t*t + tp*tp)/10.;     // exactly linear
        G_here   = 1.;     // exactly linear
        Phi_here = (t*t+tp*tp*2.)/10.;     // exactly linear
        //Phi_here = 0.9 + 2.0*t - 1.5*tp;     // exactly linear
        //Phi_here = 0.0 + 0.001*t - 0.001*tp;     // exactly linear

    std::vector<double> Vector(2);
    Vector[0] = G_here;
    Vector[1] = Phi_here;
    //std::cout << "k: " << k << ", alpha: " << alpha << ", p: " << p << " G_here  " << G_here << "  Phi_here  " << Phi_here << std::endl;
    return Vector;
}

// Generate points computes the values of both amplitude and phase functions in the requested point in (x,y) or (t,tp) plane.
// A user can chose whether a test or a real physical function is to be chosen.
std::vector<double> generatePoints(int x, int y) {
    if (test_of_phys == 1){
        return generatePoints_phys(x, y);
    }
    else if (test_of_phys == 0){
        return generatePoints_test(x, y);
    }
    else {
        std::cout << "FATAL ERROR: MODEL UNIDENTIFIED (FUNCTION GENERATEPOINTS)" << std::endl;
        exit(1);
    }
}

//****************************************************************************************************************************************
//****** END: PREPARATION OF THE FUNCTION TO BE INTEGRATED *******************************************************************************
//****************************************************************************************************************************************


//****************************************************************************************************************************************
//****** NUMERICAL SOLVER ROUTINES *******************************************************************************************************
//****************************************************************************************************************************************


// Function to compute the coefficients p, q, r in the plane equation z = px + qy + r
// once the plane is defined by the vertices of a triangle in 3D space. In our case the 3D space is constructed by the two coordinates -
// time of ionization and time of recollision - and the value of a function - the amplitude function or the phase function.
// Compute p, q, r for z(x,y) = p*x + q*y + r from three points.
// Addresses:
//   #1 Tolerant plane fit (centroid shift + C tolerance; near-vertical → constant z)
//   #5 Flat-z special case uses a tolerance instead of a hard equality test.
// ===================================================================================
// GAUSS-LEGENDRE INTEGRATOR
// ===================================================================================

// Funkcja wyznaczająca węzły (x) i wagi (w) kwadratury Gaussa-Legendre'a
void generuj_wezly_i_wagi_gaussa(int N, std::vector<double>& x, std::vector<double>& w) {
    x.resize(N);//ustalamy rozmiar węzłów/wag
    w.resize(N);//N to rząd kwadratury
    int m = (N+ 1) / 2;

    for (int i = 0; i < m; ++i) {
        double z = std::cos(M_PI * (i + 0.75) / (N + 0.5));//przybliżenie asymptotyczne Tricomieg'o- szukamy numerycznie miejsc zerowych wielomianów Legendre'a
        double z1; //poprzednie miejsce zerowe
        double p1; //wartosc wielomianu Legendre'a w z dla N
        double p2; //wartosc wielomianu Legendre'a w z dla N-1
        double p3; //wartosc wielomianu Legendre'a w z dla N-2
        double pp; //pochodna wielomianu Legendre'a w z

        do {
            p1 = 1.0; p2 = 0.0;
            for (int j = 0; j < N; ++j) {
                p3 = p2; p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
            }
            pp = N * (z * p1 - p2) / (z * z - 1.0);//obliczamy pochodna wielomianu Legendre'a
            z1 = z;
            z = z1 - p1 / pp;//obliczamy nowe znaczenie z stosujac tzw. Krok Newton'a: "szukamy dokładniejszego" miejsca zerowego z tak że z=z1-funkcja/pochodnafunkcji
        } while (std::abs(z - z1) > 1e-15);//warunek z jaka dokładnością znajdujemy miejsca zerowe

        x[i] = -z;
        x[N - 1 - i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
        w[N - 1 - i] = w[i];
    }
}

// Nowa wersja funkcji całkowania 2D oparta na metodzie Gaussa-Legendre'a
Complex integrate2D() {
    Complex integral(0.0, 0.0);

    // Generowanie n węzłów i wag w przedziale [-1, 1]
    std::vector<double> nodes, weights;
    generuj_wezly_i_wagi_gaussa(n, nodes, weights);

    // Przeskalowanie przedziału z [-1, 1] na czas [t0, tmax]
    double t_mid = 0.5 * (tmax + t0);
    double t_half = 0.5 * (tmax - t0);

    #pragma omp parallel for reduction(+:integral) collapse(2)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {

            // zmiana węzłów Gaussa na rzeczywisty czas t (rescattering) i tp (ionization)
            double t  = t_mid + t_half * nodes[i];
            double tp = t_mid + t_half * nodes[j];

            // Odcięcie obszaru poza warunkami fizycznymi- interesuje nas tp>t
            if (zero_above_diagonal == 1) {
                if (tp > t || tp > tlength) {
                    continue;
                }
            }

            // Przeliczenie na najbliższy indeks siatki dyskretnej dla tablic
            int x_idx = std::round((t - t0) / dt);
            int y_idx = std::round((tp - t0) / dt);

            // Zabezpieczenie przed wyjściem poza rozmiar tablic
            if (x_idx >= n) x_idx = n - 1;
            if (y_idx >= n) y_idx = n - 1;
            if (x_idx < 0)  x_idx = 0;
            if (y_idx < 0)  y_idx = 0;

            // Pobranie Amplitudy i Fazy dla wyznaczonego punktu
            std::vector<double> val = generatePoints(x_idx, y_idx);
            Complex f_val = val[0] * std::exp(Iu * val[1]);

            // Mnożenie przez wagę 2D (waga_i * waga_j)
            double w_2d = weights[i] * weights[j];
            integral += f_val * w_2d;
        }
    }

    // Uwzględnienie jakobianu zmiany zmiennych dziedziny-przeszliśmy na czas t i tp
    integral *= (t_half * t_half);

    std::cout << std::endl;
    return integral;
}

//****************************************************************************************************************************************
//****** END: NUMERICAL SOLVER ROUTINES *******************************************************************************************************
//****************************************************************************************************************************************


double crosssection(double p, int model) {
    //double temp = (2.0*pow(alpha, 2.0) + 4.0 * pow(p, 2.0)) / (pow(alpha, 2.0) + 4.0 * pow(p, 2.0));
    double temp;
    if (model == 1){
        temp = alpha / (pow(alpha, 2.0) + 4.0 * pow(p, 2.0));
        temp = pow(temp, 2.0);
    }
    else if (model == 2){
        temp = alpha / (pow(alpha, 2.0) + 4.0 * pow(p, 2.0));
    }
    else {
        std::cout << "FATAL ERROR: MODEL UNIDENTIFIED (FUNCTION crosssection)" << std::endl;
        exit(1);
    }
    return temp;
}

int main(int argc, char *argv[]) {
    char filename[1000];

    Complex result;
    //F0 = 0.1;

    zero_above_diagonal = 0; // 0 = genuine integration of the function; 1 = remove half of integration area, excluding the areas where the function is zero

    test_of_phys = 1; //0 = test, 1 =phys

    Ip_ion = atof(argv[1]); //ionization potential of an ion (to which the rescattering is to happen!)
    nc = atoi(argv[2]); //number of cycles of laser pulse
    omega = atof(argv[3]); //laser frequency in a.u.
    F0 = atof(argv[4]); //laser amplitude in a.u.
    model = atoi(argv[5]); //model type: 1 - 1D, 2 - Hydrogen 3D, 3 - HF 3D

    F0_nominal = F0;
    //F0 *= sqrt(2.0/3.0);    // taking into account geometry of the numerical model
    cout << F0 << endl;
    //omega = 0.06;
    //Ip = 0.39;
    Ip = 0.5;
    pi = M_PI;   // M_PI is defined in math.h
    PI=pi;
    if (model == 1){
        alpha = 1.55; // elastic scattering parameter (can be arbitrary)
        beta_par = 0.75; // initial wavefunction parameter
        epsilon = 1e-4; // 2025-06-17 A special check of dependence of the epsilon parameter has been performed for situation when sudden jump happened (p=0.27, n=20000)
        // The results have shown to be very stable for the plateau epsilon \in (1E-5,1E-3). Here we choose the value in the middle of this plateau.
    }
    else if (model ==2){
        alpha = 1.0; // elastic scattering parameter [taken as in Morishita 2007]
        beta_par = 1.0; // initial wavefunction parameter = same as nuclear charge
        epsilon = 1e-1; // like in Morishita paper (2008)
    }
    //epsilon = 1e-10;
    //epsilon = 1e-4; // 2025-06-17 A special check of dependence of the epsilon parameter has been performed for situation when sudden jump happened (p=0.27, n=20000)
    // The results have shown to be very stable for the plateau epsilon \in (1E-5,1E-3). Here we choose the value in the middle of this plateau.
    //epsilon = 1e-4; // like in Morishita paper (2008)

    period = 2 * pi / omega;
    tlength = nc * period; // length of the laser pulse



    int calculation_regime = 2;
    if (calculation_regime ==1) {
                n = int(double(nc+1)/6.0*500.0*pow((2.0*F0/omega - 0.6),2.0)); // ensures precision for momenta up to cutoff
                pmax = F0/omega*sqrt(3.17/2.0) + 0.2; //for momenta up to cutoff
    }
    else if (calculation_regime ==2) {
                n = int(double(nc+1)/6.0*500.0*pow((0.6*F0/omega + sqrt(2.0*Ip_ion) + 0.9),2.0)); // ensures precision for practically importnt momenta up to (corresponding to ionization potential of the ion onto which rescatters + 1.5 a.u.)
                pmax = sqrt(2.0*Ip_ion) + 2.5; //for practically importnt momenta up to (corresponding to ionization potential of the ion onto which rescatters + 1.5 a.u.)
    }

// n is the number of points, the number of intervals is (n-1)
n=300;

pmax=4.0;

    //tmax = 500.0;
    tmax = tlength*1.5;
    xmax = tmax;
    t0 = 0.0;
    dt = tmax/(n-1);

    //dt = xmax/n;

if (test_of_phys == 1){
    A_array.resize(n);
    F_array.resize(n);
    S1_array.resize(n);
    S2_array.resize(n);
    k_array.resize(n*n);
    check_array.resize(n*n);

    compute_A_array();
    compute_F_array();
    compute_S1_and_k_array();
    compute_S2_array();
    checking_action();
    /*
    exportVectorToFile(A_array, "A_array.txt");
    exportVectorToFile(F_array, "F_array.txt");
    exportVectorToFile(S1_array, "S1_array.txt");
    exportVectorToFile(S2_array, "S2_array.txt");
    export2DArray(k_array, n, n, "array2D.txt");
    export2DArray(check_array, n, n, "check2D.txt");
    */

}


    // Compute the integral
    //Complex result = integrate2D(x_min, x_max, y_min, y_max, n, p);
    //double upperlimit = 24/0.2*F0;
/*
    for (double p = 0.0; p <= 16.0; p += 0.01){
        if (F0 == 0.1) {
            double cs = 0;
            cs = crosssection(p * 0.55);
            std::string baseFilename = "elastic_CSs";
            // Create filename with double value
            std::ostringstream filenameStream;
            filenameStream << baseFilename << ".dat";
            std::string filename = filenameStream.str();
            writeNumbersToFile(filename, p * 0.55, cs);
        }
    }
    exit(0);
*/



   // for (double p = 0.0; p <= 16.0; p += 0.01){
    std::cout << "Calculations for field F0 = : " << F0  << std::endl;
    std::cout << "P max is : " << pmax  << std::endl;
    std::cout << "time scale is divided into n = : " << n  << std::endl;
    //for (p_current = 0.0; p_current <= pmax; p_current += 0.1){
    for (p_current = 0.0; p_current <= pmax; p_current += 0.001){
        //result = 1.0;
        std::cout << "P is : " << p_current  << std::flush;
        result = integrate2D();
        double modulus = abs(result) * abs(result);
        //modulus *= crosssection(p_current * 0.55,model);


        std::string baseFilename = "WavepacketF";

        // Create filename with double value
        std::ostringstream filenameStream;
        filenameStream << baseFilename << F0_nominal*100  << "-n=" << n << ".dat";
        std::string filename = filenameStream.str();

        writeNumbersToFile(filename, p_current * 0.55, modulus);
    }

    // Output the result
    std::cout << "The integral of f(x, y) over the region is: " << result << std::endl;

    return 0;
}




//  LOG:

//2025-10-31 Numerical solver is separated from defining a function to be integrated. The integrator is cleanedd and fixed: formula up abd down are rewritten from scratch after analytic solutions from mathematica.
//2025-11-01 exceptions for very small denominators in formulaUp and formulaDown are fixed. Now for G(x,y) on different test functions the code works correctly. space division has little impact.
// the problem persists for non-uniform G(x,y) - the results are incorrect. To be fixed.
//2025-11-02 numerical solver is now fine. Experimentally, the grid nodes number for good precision is n = 2.0/5.0 * l^2 * f'', where f'' is the second derivative (was a constant in the experiments) and l is the length of axis (time axis)
// the experiment was performed for G(x,y)=1.0 and thus does not include its influence.
// NOTE: in formulas Up and Down, for special cases (low |a| or |b| or |a-b|) the simple expressions are used (a=0 ,b=0 and (a-b)=0 correspondingly), even though these quantities may ne low, but non zero.
// In a file formulas_more_sophisticated_attempt_2025-11-01.nb the Maclaurin expansions for all those cases around small parameter are given. One can think of implementing them instead for enhancing the computational algorithm convergence speed.
// 2025-11-05 Bug found: if integrating with physocal integrated function and not only over the lower triangle but everywhere, -nan for n=1000 is revealed as a final results
// czeck k_array - indices. are they ok? Still getting -nan
