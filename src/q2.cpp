#include <iostream>
#include <iomanip>
#include <valarray>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include "q234.hpp"
#include "interp.hpp"

using namespace std;

/* Read a file with position timeseries data formatted
 * in 3 space-separated columns:
 * t x y
 * and return the numerical data in the form of a single
 * valarray of size N*3. The data is read line-by line,
 * so that the value of row i and column j is stored in
 * the (3*i + j)-th component of the valarray
 */
valarray<float> read_tracking_data(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return valarray<float>();
    }

    vector<float> data; // Use vector so we can .push_back()

    string line;
    int row = 0;
    while (getline(file, line)) {
        stringstream ss(line);
        string cell;
        int col=0;
        while (getline(ss, cell, ' ') && col < 3) {
            try {
                data.push_back(stod(cell));
            } catch (const invalid_argument& e) {
                cerr << "Error: Invalid data in file " << filename << " at row " << row << endl;
                return valarray<float>();
            }
            ++col;
        }
        ++row;
    }

    // Copy data of vector into valarray and return it
    valarray<float> valdata(data.data(), data.size());

    return valdata;
}

/* Add your functions here */

/* Transforms the tracking data into physical coordinates in metres */
void coord_xfm(valarray<float> &x_out, valarray<float> &y_out,
               valarray<float> x_data, valarray<float> y_data) {
  /* Coordinates of the center of the pitch */
  float x_center = PITCH_L / 2.0;
  float y_center = PITCH_W / 2.0;

  /* Scale and shift the normalised coordinates */
  x_out = (x_data * PITCH_L) - x_center;
  y_out = (y_data * PITCH_W) - y_center;
}

/* Compute Centered finite differences:
     * vx[i] = (x[i+1] - x[i-1])/(2*dt)
     * vy[i] = (y[i+1] - y[i-1])/(2*dt) */
void comp_speed(valarray<float> &vx, valarray<float> &vy,
    const valarray<float> &x_data, const valarray<float> &y_data,
    float dt) {
    const size_t N = x_data.size() - 2;
    vx.resize(N);
    vy.resize(N);

    vx = (x_data[slice(2, N, 1)] - x_data[slice(0, N, 1)]) / (2.0f * dt);
    vy = (y_data[slice(2, N, 1)] - y_data[slice(0, N, 1)]) / (2.0f * dt);
}

/* Compute the acceleration using second-order accurate centered difference for second order derivative:
* ax[i] = (x[i+1] - 2 * x[i] + x[i-1])/(dt^2)
* ay[i] = (y[i+1] - 2 * y[i} + y[i-1])/(dt^2) */
void comp_acc(valarray<float> &ax, valarray<float> &ay,
            const valarray<float> &x_data, const valarray<float> &y_data, float dt) {
    /* We cannot compute centered finite differences at the endpoints */
    const size_t N = x_data.size() - 2;
    ax.resize(N);
    ay.resize(N);

    ax = (x_data[slice(2, N, 1)]
            - (2.0f * x_data[slice(1, N, 1)])
            + x_data[slice(0, N, 1)]) / (dt*dt);

    ay = (y_data[slice(2, N, 1)]
            - (2.0f * y_data[slice(1, N, 1)])
            + y_data[slice(0, N, 1)]) / (dt*dt);
}

/* Compute magnitude: mag = sqrt(vx^2 + vy^2) for each pair of components */
void compute_mag(const valarray<float> &x, const valarray<float> &y,
    valarray<float>& mag) {
    mag.resize(x.size());
    mag = sqrt((x*x) + (y*y));
}

int main(int argc, char *argv[]) {

    valarray<float> data = read_tracking_data("tracking_data.dat");
    size_t Ndata = data.size()/3;

    /* Initialise new individual arrays by taking slices of
     * data (with stride=3) to extract the three columns  */
    valarray<float> t_varr(data[slice(0, Ndata, 3)]);
    valarray<float> x_varr(data[slice(1, Ndata, 3)]);
    valarray<float> y_varr(data[slice(2, Ndata, 3)]);

    /* Continue the main() here  */

    /* Transform the tracking data to physical coordinates */
    valarray<float> x_physical(Ndata), y_physical(Ndata);
    coord_xfm(x_physical, y_physical, x_varr, y_varr);

    /*                             SPEED                                    */
    /* Initialise new arrays for the speed components of the player
     * Compute speed components vx and vy*/
    float dt = 0.2;
    valarray<float> vx, vy, v_mag;
    comp_speed(vx, vy, x_physical, y_physical, dt);
    /* Compute the speed magnitude and store it in v_mag */
    compute_mag(vx, vy, v_mag);

    /*                         ACCELERATION                                */
    /* Initialise new arrays for the acceleration components of the player
     * and compute acceleration components ax and ay*/
    valarray<float> ax, ay, a_mag;
    comp_acc(ax, ay, x_physical, y_physical, dt);
    /* Compute the acceleration magnitude and store it in a_mag */
    compute_mag(ax, ay, a_mag);

    /* Write out file */
    ofstream speed_data("speed_data.dat");
    speed_data << fixed << setprecision(6);
    for (size_t i = 0; i < v_mag.size(); ++i) {
        speed_data << t_varr[i+1] << "\t" << v_mag[i] << endl;
    }

    /* Find maximum speed */
    cout << "Maximum speed = " << v_mag.max() << " m/s" << endl;

    return 0;
}
