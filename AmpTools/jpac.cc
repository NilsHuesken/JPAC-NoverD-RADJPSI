#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include "TMath.h"

#include "IUAmpTools/Kinematics.h"
#include "PROJECTAmp/jpac.h"

jpac::jpac( const vector< string >& args ) :
UserAmplitude< jpac >(args)
{
  m_p1 = atoi(args.at(0).c_str());
  m_p2 = atoi(args.at(1).c_str());
  m_J = atoi(args.at(2).c_str());
  m_idx = atoi(args.at(3).c_str());
  m_phsp = args.at(4);

  v_m.push_back(0.1349766);
  v_m.push_back(0.4976);
  v_m.push_back(0.762);

  string radjpsi = std::getenv("PROJECT");
  string filename = radjpsi + "/PROJECTExe/PHSPFUNCS/rhoN-" + m_phsp + "_bootstrap.txt";
  ifstream in(filename.c_str());
  assert(in.is_open());
  string line;
  while(in.is_open() && getline(in,line)){
    string s[13];
    istringstream ss(line);
    ss >> s[0] >> s[1] >> s[2] >> s[3] >> s[4] >> s[5] >> s[6] >> s[7] >> s[8] >> s[9] >> s[10] >> s[11] >> s[12];
    double sqrts = atof(s[0].c_str());
    complex<double> ch1, ch2, ch3;
    if(m_J==0){
      ch1 = complex<double>(atof(s[1].c_str()),atof(s[2].c_str()));
      ch2 = complex<double>(atof(s[3].c_str()),atof(s[4].c_str()));
      ch3 = complex<double>(atof(s[5].c_str()),atof(s[6].c_str()));
    }
    else{
      ch1 = complex<double>(atof(s[7].c_str()),atof(s[8].c_str()));
      ch2 = complex<double>(atof(s[9].c_str()),atof(s[10].c_str()));
      ch3 = complex<double>(atof(s[11].c_str()),atof(s[12].c_str()));
    }
    v_int[0].push_back(make_pair(sqrts,ch1));
    v_int[1].push_back(make_pair(sqrts,ch2));
    v_int[2].push_back(make_pair(sqrts,ch3));
  }

  if(m_phsp.compare("inputcddcm11newc3")==0){ // 1
      if(m_J==0){ 
        m_model_K = 1;
        m_model_n = 2;
      }
      else{
        m_model_K = 0;
        m_model_n = 2;
      }
  }
  else if(m_phsp.compare("inputcddcm3newc3")==0){ // 2
      if(m_J==0){ 
        m_model_K = 1;
        m_model_n = 1;
      }
      else{
        m_model_K = 0;
        m_model_n = 1;
      }
  }
  else if(m_phsp.compare("inputcddcm5newc3")==0){ // 3
      if(m_J==0){ 
        m_model_K = 1;
        m_model_n = 2;
      }
      else{
        m_model_K = 0;
        m_model_n = 2;
      }
  }
  else if(m_phsp.compare("inputcddcm6newc3")==0){ // 4
      if(m_J==0){ 
        m_model_K = 1;
        m_model_n = 2;
      }
      else{
        m_model_K = 0;
        m_model_n = 2;
      }
  }
  else if(m_phsp.compare("inputcddcm7newc3")==0){ // 5
      if(m_J==0){ 
        m_model_K = 1;
        m_model_n = 0;
      }
      else{
        m_model_K = 0;
        m_model_n = 0;
      }
  }
  else if(m_phsp.compare("inputcddcm9newc3")==0){ // 6
      if(m_J==0){ 
        m_model_K = 1;
        m_model_n = 1;
      }
      else{
        m_model_K = 0;
        m_model_n = 1;
      }
  }
  else if(m_phsp.compare("inputcddcmnewc3")==0){ // 7
      if(m_J==0){ 
        m_model_K = 1;
        m_model_n = 0;
      }
      else{
        m_model_K = 0;
        m_model_n = 0;
      }
  }
  else if(m_phsp.compare("inputstartcm11new2c3")==0){ // 8
      if(m_J==0){ 
        m_model_K = 0;
        m_model_n = 2;
      }
      else{
        m_model_K = 0;
        m_model_n = 2;
      }
  }
  else if(m_phsp.compare("inputstartcm133c3")==0){ // 9
      if(m_J==0){ 
        m_model_K = 0;
        m_model_n = 2;
      }
      else{
        m_model_K = 0;
        m_model_n = 2;
      }
  }
  else if(m_phsp.compare("inputstartcm15new2c3")==0){ // 10
      if(m_J==0){ 
        m_model_K = 0;
        m_model_n = 2;
      }
      else{
        m_model_K = 0;
        m_model_n = 2;
      }
  }
  else if(m_phsp.compare("inputstartcm3c3")==0){ // 11
      if(m_J==0){ 
        m_model_K = 0;
        m_model_n = 1;
      }
      else{
        m_model_K = 0;
        m_model_n = 1;
      }
  }
  else if(m_phsp.compare("inputstartcm5newc3")==0){ // 12
      if(m_J==0){ 
        m_model_K = 0;
        m_model_n = 1;
      }
      else{
        m_model_K = 0;
        m_model_n = 1;
      }
  }
  else if(m_phsp.compare("inputstartcm6new2c3")==0){ // 13
      if(m_J==0){ 
        m_model_K = 0;
        m_model_n = 0;
      }
      else{
        m_model_K = 0;
        m_model_n = 0;
      }
  }
  else if(m_phsp.compare("inputstartcmnew2c3")==0){ // 14
      if(m_J==0){ 
        m_model_K = 0;
        m_model_n = 0;
      }
      else{
        m_model_K = 0;
        m_model_n = 0;
      }
  }
  




  m_a[0][0] = AmpParameter(args.at(4+1)); //-0.170; // pipi
  m_a[0][1] = AmpParameter(args.at(5+1)); //-3.471;
  m_a[0][2] = AmpParameter(args.at(6+1)); //1.087;
  m_a[0][3] = AmpParameter(args.at(7+1)); //-0.198;

  m_a[1][0] = AmpParameter(args.at(8+1)); //4.705; // KK
  m_a[1][1] = AmpParameter(args.at(9+1)); //6.541;
  m_a[1][2] = AmpParameter(args.at(10+1)); //-2.844;
  m_a[1][3] = AmpParameter(args.at(11+1)); //-0.893;

  m_a[2][0] = AmpParameter(args.at(12+1)); //0; // rho rho
  m_a[2][1] = AmpParameter(args.at(13+1)); //0;
  m_a[2][2] = AmpParameter(args.at(14+1)); //0;
  m_a[2][3] = AmpParameter(args.at(15+1)); //0;


  m_g[0][0] = AmpParameter(args.at(16+1)); //-2.127; // pipi
  m_g[1][0] = AmpParameter(args.at(17+1)); //0.779;
  m_g[2][0] = AmpParameter(args.at(18+1)); //-11.473;

  m_g[0][1] = AmpParameter(args.at(19+1)); //-1.110; // KK
  m_g[1][1] = AmpParameter(args.at(20+1)); //3.257;
  m_g[2][1] = AmpParameter(args.at(21+1)); //-6.087;

  m_g[0][2] = AmpParameter(args.at(22+1)); //0.371; // rho rho
  m_g[1][2] = AmpParameter(args.at(23+1)); //-0.236;
  m_g[2][2] = AmpParameter(args.at(24+1)); //-3.125;



  m_m0[0] = AmpParameter(args.at(25+1)); //sqrt(2.140);
  m_m0[1] = AmpParameter(args.at(26+1)); //sqrt(2.468);
  m_m0[2] = AmpParameter(args.at(27+1)); //sqrt(4.788);


  m_c[0][0] = AmpParameter(args.at(28+1)); //-8.984;
  m_c[0][1] = AmpParameter(args.at(29+1)); //4.732;
  m_c[0][2] = AmpParameter(args.at(30+1)); //-25.612;
  m_c[1][1] = AmpParameter(args.at(31+1)); //21.301;
  m_c[1][2] = AmpParameter(args.at(32+1)); //-6.836;
  m_c[2][2] = AmpParameter(args.at(33+1)); //0.489;


  m_d[0][0] = AmpParameter(args.at(34+1)); //-25.360;
  m_d[0][1] = AmpParameter(args.at(35+1)); //-14.444;
  m_d[0][2] = AmpParameter(args.at(36+1)); //0;
  m_d[1][1] = AmpParameter(args.at(37+1)); //-8.864;
  m_d[1][2] = AmpParameter(args.at(38+1)); //0;
  m_d[2][2] = AmpParameter(args.at(39+1)); //0;

  for(int i=0;i<3;i++){
    for(int j=0;j<4;j++){
      registerParameter(m_a[i][j]);
    }
  }
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      registerParameter(m_g[i][j]);
    }
  }
  for(int j=0;j<3;j++){
    registerParameter(m_m0[j]);
  }
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      if(j>=i){
        registerParameter(m_c[i][j]);
      }
    }
  }
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      if(j>=i){
        registerParameter(m_d[i][j]);
      }
    }
  }

  if(args.size()>41)
    m_L = atoi(args.at(41).c_str());
}

void
jpac::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  TLorentzVector P1, P2;
  P1.SetPxPyPzE(pKin[m_p1][1], pKin[m_p1][2], pKin[m_p1][3], pKin[m_p1][0]);
  P2.SetPxPyPzE(pKin[m_p2][1], pKin[m_p2][2], pKin[m_p2][3], pKin[m_p2][0]);

  double sqrts = (P1+P2).M();
  double s = (P1+P2).M2();
  double s0 = 1;
  double smin = pow(1,2);
  double smax = pow(2.5,2);

  userVars[ks] = s;

  double Egam = (pow(3.0969,2)-s)/(2*sqrts);
  VectorXd p = VectorXd::Zero(3);
  for(int i=0;i<3;i++)
    p(i) = sqrt(s-4*pow(v_m.at(i),2))/2.;

  userVars[kEgam] = Egam;
  userVars[kp0] = p(0);
  userVars[kp1] = p(1);
  userVars[kp2] = p(2);

  VectorXd n_vec = VectorXd::Zero(3);
  
  double omega = 0;
  if(m_model_n==0)
    omega = omega_pole(s,s0);
  else if(m_model_n==1)
    omega = omega_scaled(s,smin,smax);
  else if(m_model_n==2)
    omega = omega_polescaled(s,smin,smax,s0);
  
  userVars[kCheby0] = Cheby(omega,0);
  userVars[kCheby1] = Cheby(omega,1);
  userVars[kCheby2] = Cheby(omega,2);
  userVars[kCheby3] = Cheby(omega,3); 

  MatrixXcd integral = MatrixXcd::Zero(3,3);
  integral(0,0) = findInt(sqrts,0,v_int);
  integral(1,1) = findInt(sqrts,1,v_int);
  integral(2,2) = findInt(sqrts,2,v_int);

  userVars[kInt0re] = integral(0,0).real();
  userVars[kInt0im] = integral(0,0).imag();
  userVars[kInt1re] = integral(1,1).real();
  userVars[kInt1im] = integral(1,1).imag();
  userVars[kInt2re] = integral(2,2).real();
  userVars[kInt2im] = integral(2,2).imag();
}


complex< GDouble >
jpac::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

  double s = userVars[ks];
  double sqrts = sqrt(s);
  double Egam = userVars[kEgam];
  VectorXd p = VectorXd::Zero(3);
  p(0) = userVars[kp0];
  p(1) = userVars[kp1];
  p(2) = userVars[kp2];
  VectorXd n_vec = VectorXd::Zero(3);
  n_vec(0) += (m_a[0][0]*userVars[kCheby0]+m_a[0][1]*userVars[kCheby1]+m_a[0][2]*userVars[kCheby2]+m_a[0][3]*userVars[kCheby3]);
  n_vec(1) += (m_a[1][0]*userVars[kCheby0]+m_a[1][1]*userVars[kCheby1]+m_a[1][2]*userVars[kCheby2]+m_a[1][3]*userVars[kCheby3]);
  n_vec(2) += (m_a[2][0]*userVars[kCheby0]+m_a[2][1]*userVars[kCheby1]+m_a[2][2]*userVars[kCheby2]+m_a[2][3]*userVars[kCheby3]);

  MatrixXcd integral = MatrixXcd::Zero(3,3);
  integral(0,0) = complex<double>(userVars[kInt0re],userVars[kInt0im]);
  integral(1,1) = complex<double>(userVars[kInt1re],userVars[kInt1im]);
  integral(2,2) = complex<double>(userVars[kInt2re],userVars[kInt2im]);
  int n_res = 3;
  MatrixXd K_inv = MatrixXd::Zero(3,3);
  if(m_model_K==0){
    MatrixXd K = MatrixXd::Zero(3,3);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        if(i>j){
          K(i,j) = K(j,i);
          continue;
        }
        K(i,j) = m_c[i][j] + m_d[i][j]*s;
        for(int r=0;r<3;r++)
          K(i,j) += m_g[r][i]*m_g[r][j]/(pow(m_m0[r],1)-s);
      }
    }
    K_inv = K.inverse();
  }
  else{
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        if(i>j){
          K_inv(i,j) = K_inv(j,i);
        }
        else{
          K_inv(i,j) = pow(m_c[i][j],2) - pow(m_d[i][j],2)*s; // these parameters are squared in reality, but not in the paper!
          for(int r=0;r<3;r++){
            K_inv(i,j) -= m_g[r][i]*m_g[r][j]/(pow(m_m0[r],1)-s);
          }
        }
      }
    }
  }
  MatrixXcd K_invc = K_inv.cast<std::complex<double>>();
  MatrixXcd D = K_invc - integral;
  MatrixXcd D_inv = D.inverse();
  MatrixXd _ISO = MatrixXd::Zero(3,3);
  _ISO(0,0) = sqrt(3.);
  _ISO(1,1) = sqrt(4.);
  _ISO(2,2) = sqrt(3.);
  MatrixXcd D_ISO = D_inv * _ISO;
  

  complex<double> amp = complex<double>(0.,0.);
  for(int k=0;k<3;k++){
    amp += n_vec(k)*D_ISO(k,m_idx); 
  }
  amp *= sqrt(p(m_idx))*pow(Egam,m_L)*pow(p(m_idx),m_J);

  if(TMath::IsNaN(amp.real()) || TMath::IsNaN(amp.imag())){
    cout << __FILE__ << " " << sqrts << " " << amp << " " << p(m_idx) << " " << m_idx << endl;
    cout << "n " << n_vec << endl;
    cout << "D " << D_ISO << endl;
    for(int i=0;i<3;i++){
      cout << "m" << i << " " << m_m0[i] << endl;
      for(int j=0;j<3;j++){
        cout << "g" << i << j << " " << m_g[i][j] << endl;
        if(i>j)
          continue;
        cout << "c" << i << j << " " << m_c[i][j] << endl;
        cout << "d" << i << j << " " << m_d[i][j] << endl;
        
      }
    }
    assert(0);
  }
  
  return amp;
}



#ifdef GPU_ACCELERATION
void
jpac::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
  
  jpac_exec( dimGrid,  dimBlock, GPU_AMP_ARGS );

}
#endif //GPU_ACCELERATION
