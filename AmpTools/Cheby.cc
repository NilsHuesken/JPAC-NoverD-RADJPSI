
#include "./Cheby.h"

#include <math.h>

double Cheby(double x, int n){
  if(n==0){
    return 1;
  }
  else if(n==1){
    return x;
  }
  else{
    return 2*x*Cheby(x,n-1)-Cheby(x,n-2);
  }
}

double omega_pole(double s, double s0){
  return s/(s+s0);
}

double omega_scaled(double s, double smin, double smax){
  return 2*(s-smin)/(smax-smin)-1;
}

double omega_polescaled(double s, double smin, double smax, double s0){
  return 2*(omega_pole(s,s0)-omega_pole(smin,s0))/(omega_pole(smax,s0)-omega_pole(smin,s0))-1;
}

complex<double> findInt(double sqrts, int chanID, map<int,vector<pair<double,complex<double>>>> v_int){

  vector<pair<double,complex<double>>> vec = v_int[chanID];
  assert(vec.size());

  complex<double> returnMe;

  for(unsigned int i=1;i<vec.size();i++){
    if(vec.at(i).first > sqrts){
      double this_sqrts = vec.at(i).first;
      double last_sqrts = vec.at(i-1).first;
      double this_int_re = vec.at(i).second.real();
      double this_int_im = vec.at(i).second.imag();
      double last_int_re = vec.at(i-1).second.real();
      double last_int_im = vec.at(i-1).second.imag();
      double m_re = (this_int_re - last_int_re)/(this_sqrts - last_sqrts);
      double b_re = this_int_re - m_re * this_sqrts;
      double m_im = (this_int_im - last_int_im)/(this_sqrts - last_sqrts);
      double b_im = this_int_im - m_im * this_sqrts;
      returnMe = complex<double>(m_re*sqrts+b_re,m_im*sqrts+b_im);
      break;
    }
  }
  return returnMe;
}
