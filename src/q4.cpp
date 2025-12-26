#include <iostream>
#include <valarray>
#include <cmath>
#include "q234.hpp"

using namespace std;

/* Runge-Kutta 4th order method */
valarray<float> rk4(float t0, valarray<float> y0, float h, int N) {
  float t = t0;
  valarray<float> y = y0;
  valarray<float> k1, k2, k3, k4;

  // Evolution loop for 4th order Runge-Kutta
  for (int i=0; i<N; i++) {
    k1 = h * rhs(t, y);
    k2 = h * rhs(t + h / 2.0, y + k1 / 2.0);
    k3 = h * rhs(t + h / 2.0, y + k2 / 2.0);
    k4 = h * rhs(t + h, y + k3);

    // Update y using the weighted average of slopes
    y += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

    // Move to the next step
    t += h;
  }

  return y;
}

/* Returns the derivative */
valarray<float> rhs(float t, const valarray<float> yvec){
  valarray<float> dydt(6);
  dydt[0] = yvec[3]; /* x' = vx */
  dydt[1] = yvec[4]; /* y' = vy */
  dydt[2] = yvec[5]; /* z' = vz */

  /* If drag_on is true then define d = C_d*rho*A/2m if false = 0*/
  float d = drag_on ? 0.5 * C_DRAG * RHO_AIR * M_PI * R_BALL * R_BALL / M_BALL : 0;

  /* If magnus_on is true then define magnus, if false = 0 */
  float magnus = magnus_on ? S_MAGN * omega[2] / M_BALL : 0;

  /* Determine |v| = sqrt(vx^2 + vy^2 + vz^2) */
  float vmag = sqrt(yvec[3]*yvec[3] + yvec[4]*yvec[4] + yvec[5]*yvec[5]);

  /* Acceleration components */
  dydt[3] = -d * vmag * yvec[3] - magnus * yvec[4];
  dydt[4] = -d * vmag * yvec[4] + magnus * yvec[3];
  dydt[5] = -A_GRAV - d * vmag * yvec[5];

  return dydt;
}

int main(int argc, char *argv[]){

  /* Initial time and position */
  float t0 = 0.0;
  float x0 = PITCH_L * 0.5 - 20;/* 20m in front of the goal line */
  float y0 = 3.0;               /* 3m left of goal centre */
  float z0 = R_BALL;            /* touching the ground */

  /* Initial velocity and angle */
  float v0 = 25.0;                        /* Default Launch velocity */
  if (argc > 1){ v0 = stof(argv[1]); } /* If the value is given override */
  float theta = 20*M_PI/180;              /* Launch agle in radians */

  /* Initial velocity components */
  float vx0 = v0 * cos(theta); /* x component */
  float vy0 = 0.0;             /* No sideways component */
  float vz0 = v0 * sin(theta); /* z component */

  /* Initialise the array Y with values at t=0 */
  valarray<float> Y = {x0, y0, z0, vx0, vy0, vz0};

  cout << fixed << setprecision(6);
  cout << t0 << ' '
  << Y[0] << ' ' << Y[1] << ' ' << Y[2] << ' '
  << Y[3] << ' ' << Y[4] << ' ' << Y[5] << endl;

  float dt = 0.01; /* Time step */

  /* If the ball's centre hasn't yet reached the goal line
   * and the ball's centre is still above the ground */
  while ( Y[0] < 0.5*PITCH_L && Y[2] > R_BALL) {
    Y = rk4(t0, Y, dt, 1); /* Update the state of Y using Runge Kutta */
    t0 += dt;                      /* Increase time by dt */

    /* Each row of the output has t x y z vx vy vz */
    cout << t0 << ' '
    << Y[0] << ' ' << Y[1] << ' ' << Y[2] << ' '
    << Y[3] << ' ' << Y[4] << ' ' << Y[5] << endl;
  }
  if (Y[0] >= 0.5*PITCH_L) {     /* If the ball reaches the goal line */
    float top_z = Y[2] + R_BALL; /* Distance from the top of the ball to the ground =
                                  *z position of the centre of the ball + radius */
    if ( top_z < GOAL_H )
      cerr << "GOOD SHOT!" << endl;
    else
      cerr << "The ball went over the horizontal bar" << endl;
  }
  else { /* The ball has not reached the goal line */
    cerr << "WEAK SHOT!" << endl;
  }
  return 0;


}
