#include <math.h>
#include <iostream>
#include <complex>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <stdlib.h>
#include <utility>
#include <algorithm>
#include <cassert>

using namespace std;
using std::complex;

double Cheby(double x, int n);
double omega_pole(double s, double s0);
double omega_scaled(double s, double smin, double smax);
double omega_polescaled(double s, double smin, double smax, double s0);
complex<double> findInt(double sqrts, int chanID, map<int,vector<pair<double,complex<double>>>> v_int);

