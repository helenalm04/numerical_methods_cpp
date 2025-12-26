#include <iostream>
#include <valarray>
#include <fstream> // so we can use ifstream
#include <vector>
#include <cassert>
#include <cstdlib> // so we can use exit()
#include "q234.hpp"
#include "interp.hpp"
using namespace std;

/* Calculates array of trapezoid coefficients c_i,
 * given N segments of equal step-size h. Single-step
 * corresponds to N=1, h = b-a */
valarray<float> trapz_coeffs(size_t N, float h) {
    valarray<float> ci(1.0, N+1);
    ci[0] = 0.5;
    ci[N] = 0.5;

    ci *= h;
    return ci;
}

/* Calculates vector of simpson coefficients c_i,
 * given N segments of equal step-size h. Single-step
 * corresponds to N=2, h = (b-a)/2 */
valarray<float> simpson_coeffs(size_t N, float h) {
    if (N%2) {
        cout << "ERROR: Need even number of segments for composite Simpson integration!" << endl;
        exit(-1);
    }
    valarray<float> ci(2.0/3.0, N+1);
    for (int i = 0; i < ci.size(); i++) {
        ci[i] += 2.0/3.0 * (i % 2);
    }
    ci[0] = ci[N] = 1.0/3.0;
    ci *= h;
    return ci;
}

/* Calculates vector of simpson's 3/8 coefficients c_i,
 * given N segments of equal step-size h.*/
valarray<float> simpson_four_coeffs(size_t N, float h) {
    valarray<float> ci(0.0, N+1);
    for (size_t i=0; i+3 <= N; i += 3) {
        ci[i] += 1.0;
        ci[i+1] += 3.0;
        ci[i+2] += 3.0;
        ci[i+3] += 1.0;
    }
    ci *= 3.0 * h / 8.0;
    /* Simpson's rule can be used for the remaining subintervals
     * without changing the order of the error term */
    size_t left = N%3;
    if (left != 0) {
        cerr << "Simpson's 3/8 rule: " << left
        << " leftover segment" << (left == 1 ? "" : "s")
        << " integrated with ";
        size_t s = N - left;
        if (left == 1) {
            cerr << "trapezoid rule." << endl;
            ci[s] += h / 2.0;
            ci[s+1] += h / 2.0;
        } else {
            cerr << "Simpson's 1/3 rule." << endl;
            /* two segments */
            ci[s]     += h / 3.0;   // first endpoint
            ci[s + 1] += 4.0 * h / 3.0;   // midpoint
            ci[s + 2] += h / 3.0;   // last endpoint
        }
    }
    return ci;
}

/* We define a function that takes care of numerical integration in 1D,
 * given the array with values f(x_i), the length of the intevals and
 * the Newton-Cotes formula to use, based on its number of points,
 * e.g. 2: trapezoid, 3: Simpson, etc. */
float nintegrate1D(valarray<float> f, float h, int NCpoints) {
    size_t N = f.size() - 1; /* The number of segments */

    /* Compute array of coefficients c_i based on the chosen rule */
    valarray<float> ci(N+1);

    switch (NCpoints) {
        case 2:  // trapezoid rule
            ci = trapz_coeffs(N, h);
        break;
        case 3:  // Simpson's rule
            ci = simpson_coeffs(N, h);
        break;
        case 4: // Simpson's 3/8 rule
            ci = simpson_four_coeffs(N, h);
        break;
        default:
            cout << "Newton-Cotes " << NCpoints << "-point formula not available. Exiting..." << endl;
        return -1;
    }

    /* perform discrete sum \sum_{i=0}^{N} c_i f_i */
    return (ci*f).sum();
}

int main() {
    valarray <float> v_data{0.0, 3.0, 5.0, 8.0};
    valarray <float> P_data{100.0, 700.0, 1100.0, 2000.0};

    /* Evaluate P_A at v=5.8 */
    float v_eval = 5.8;
    float La = Lagrange_N(v_data, P_data, v_eval); // defined in interp.hpp
    cout  << "Evaluating P_A(v) when v=5.8 we get: P(" << v_eval << ") = " <<  La << endl;

    /* Read speed data from speed_data.dat */
    vector<float> t_vec, v_vec; // Use vector so we can .push_back()
    ifstream file("speed_data.dat");
    if(!file) {
        cerr << "Error: Unable to open file." << endl;
        return 1;
    }
    float t,v;
    while (file>>t>>v) { /* Keep reading until the end of file */
        t_vec.push_back(t);
        v_vec.push_back(v);
    }
    size_t Ndata = v_vec.size(); /* Number of points */

    /* Convert vector into valarray so we can use it in interp_coeffs */
    valarray<float> t_varr(t_vec.data(), Ndata); /* Time data */
    valarray<float> v_varr(v_vec.data(), Ndata); /* Velocity data  */

    /* Compute coefficients */
    valarray<float> coeffs = interp_coeffs(v_data, P_data); // defined in interp.hpp

    /* Copy coefficients and speed data into its own vector so we can use them in poly_eval */
    vector<float> coeff_vec(&coeffs[0], &coeffs[0] + coeffs.size());
    vector<float> vdata_vec (&v_data[0],  &v_data[0] + v_data.size());

    /* Build P_i = P_A(v_i) */
    valarray <float> p_i(Ndata);
    for (size_t i = 0; i < Ndata; ++i) {
        p_i[i] = poly_eval(coeff_vec, vdata_vec, v_varr[i]); // defined in interp.hpp
    }
    /* Time step */
    float dt = t_varr[1] - t_varr[0];

    /* Use the trapezium composite rule to integrate P vs t to get energy in joules */
    float E_trapz = nintegrate1D(p_i, dt, 2);
    cout << "Energy (using trapezium rule) = " << E_trapz << " J" << endl;

    /* Use the Simpson's 3/8 to integrate P vs t to get energy in joules */
    float E_sim38  = nintegrate1D(p_i, dt, 4);
    cout << "Energy (using Simpson 3/8) = " << E_sim38 << " J" << endl;

    return 0;
 }


