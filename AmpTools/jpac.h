#if !defined(JPAC)
#define JPAC

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"
#include "TLorentzVector.h"
#include "PROJECTAmp/Cheby.h"
#include "../../Eigen/Dense"
#include "TString.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>
#include <algorithm>

#ifdef GPU_ACCELERATION

void jpac_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;
using namespace Eigen;

class Kinematics;

class jpac : public UserAmplitude< jpac >{

public:
  
  jpac() : UserAmplitude< jpac >() { }

  jpac( const vector< string >& args );

  ~jpac(){}

  string name() const { return "jpac"; }

  complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
  

  enum UserVars { ks=0, kEgam, kp0, kp1, kp2, kCheby0, kCheby1, kCheby2, kCheby3, kInt0re, kInt0im, kInt1re, kInt1im, kInt2re, kInt2im, kNumUserVars };
  unsigned int numUserVars() const { return kNumUserVars; }
  
  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;
  
  bool needsUserVarsOnly() const { return true; }

  //void updatePar(const AmpParameter &par);
  // **  end of optional lines **
  
#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

  bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	int m_p1, m_p2;
  vector<double> v_m;
  AmpParameter m_a[5][5];
  AmpParameter m_g[5][5];
  AmpParameter m_m0[5];
  AmpParameter m_c[5][5], m_d[5][5];

  int m_model_n = 1;
  int m_model_K = 1;
  int m_model_INT = 1;

  double m_sL = 0;
  double m_alpha = 0;
  int m_J = 0;
  int m_L = 1;
  int m_idx;
  string m_phsp;
  bool recalc_D;

  map<int,vector<pair<double,complex<double>>>> v_int;
};

#endif
