#ifndef INTERP_H_
#define INTERP_H_

#include <iostream>
#include <cassert>
#include <valarray>
#include <vector>


using namespace std;

int locate(valarray<float>& xi, float x, bool uniform=true);

/* Given an array xi of n grid points, calculate
 * the k-th Lagrange polynomial of degree n-1 that
 * has roots at each grid point except for xi[k]
 * and that is normalised so that its value on
 * xi[k] is 1.0. */
float Lagrange_Nk(int k, valarray<float>& xi, float x) {

    size_t n = xi.size();
    float prod = 1.0;

    /* For loop that calculates the product
     * of Eq.(7) from the Lecture Notes */
    for (size_t i=0;i<n;i++) {
        if (i==k) continue;
        prod *= (x - xi[i])/(xi[k] - xi[i]);
    }

    return prod;
}

/* Given an array xi of n grid points, and an array
 * yi of n values, calculate the Lagrange interpolant
 * that passes from all data points (xi[i],yi[i]) and
 * return its value at a given point x in [xi[0],xi[n-1]]. */
float Lagrange_N(valarray<float>& xi, valarray<float>& yi, float x) {

    size_t n = xi.size();
    assert(yi.size() == n);

    int k;
    float sum = 0.0;

    /* Perform the linear superposition of all L_Nk(x) */
    for (k=0;k<n;k++) {
        sum += yi[k] * Lagrange_Nk(k, xi, x);
    }

    return sum;
}

/* Alternative way of performing interpolation */

/* Compute the coefficients of the interpolating polynomial
 * for the data (x_i, y_i), using Newton's formula for
 * divided differences. */
valarray<float> interp_coeffs(const valarray<float>& xi, const valarray<float>& yi) {
    int n = xi.size();
    vector < vector<float> > divided_differences(n, vector<float>(n, 0.0));
    valarray<float> coeffs(n);

    // Initialize the first column with yi values
    for (int i = 0; i < n; i++) {
        divided_differences[i][0] = yi[i];
    }

    // Compute divided differences table
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            divided_differences[i][j] =
                (divided_differences[i + 1][j - 1] - divided_differences[i][j - 1]) / (xi[i + j] - xi[i]);
        }
    }

    /* Extract coefficients from divided differences */
    for (int i = 0; i < n; i++) {
        coeffs[i] = divided_differences[0][i];
    }

    return coeffs;
}

/* Function that evaluates the interpolating polynomial at a given x
 * given its coefficients that interp_coeffs() computed for the grid xi */
float poly_eval(const vector<float>& coeffs, const vector<float>& xi, float x) {
    int n = coeffs.size();

    /* Evaluate interpolating polynomial using nested form.
     * Start with higher order coefficient (see e.g. Numerical
     * Recipes 3rd Ed. Sec.3.2) */
    float result = coeffs[n - 1];

    for (int i = n - 2; i >= 0; i--) {
        result = result * (x - xi[i]) + coeffs[i];
    }

    return result;
}

#endif // INTERP_H_
