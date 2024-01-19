# JPAC-NoverD-RADJPSI

This code reproduces the results shown in *A. Rodas et al. (JPAC Collaboration), Eur.Phys.J.C 82 (2022) 1, 80 (DOI: 10.1140/epjc/s10052-022-10014-8, arXiv:2110.00027*), using the parameter values given in the Appendix.
The dispersion integrals in eq. (4) for the different choices of rhoN have been provided in the form of lookup-tables by Alessandro Pilloni and Arkaitz Rodas.

For matrix calculations, the code uses the *Eigen library* (see https://eigen.tuxfamily.org) - a copy of Eigen is **not** included in this repository. To run the code, you need to download a copy of Eigen and adjust the include-path for Eigen/Dense.

ROOT is used to plot the squared-magnitude of the amplitude compared to binned BESIII data for the pipi- (*Phys.Rev.D 92 (2015) 5, 052003*) and KK-channels (*Phys.Rev.D 98 (2018) 7, 072003*).

If useful, please feel free to try and use this amplitude as a description of the pipi or KK S- and D-wave between 1 and 2.5 GeV in your work, citing the original publication *Eur.Phys.J.C 82 (2022) 1, 80* as well as https://eigen.tuxfamily.org.
You can also refer to doi.org/XXXX/zenodo.XXXX for this repository.

### AmpTools integration
In the folder *AmpTools*, there is also an AmpTools version of this amplitude (see https://github.com/mashephe/AmpTools and doi.org/10.5281/zenodo.5039377 for details on AmpTools). To use it, copy it into your AmpTools project under *PROJECT/PROJECTLib/PROJECTAmp*. You then need to change the word *PROJECT* to your project-name in *jpac.cc* and *jpac.h*, make sure you place the lookup-tables in *PROJECT/PROJECTExe/PHSPFUNCS* and again adjust the include-path for Eigen/Dense.
In an AmpTools config-file, the amplitude can then be used like this:
```
amplitude reaction::sum::amp jpac i1 i2 L f lookup-table parameterlist
```
where i1 and i2 are the indices of the two pions/kaons, L=0,2 switches between S- and D-wave, f is the final-state index (0: pipi, 1: KK), lookup-table is the string between rhoN- and _bootstrap.txt in PHSPFuncs/ and parameterlist is a proxy for
```
define parameterlist [a00] [a01] [a02] [a03] [a10] [a11] [a12] [a13] [a20] [a21] [a22] [a23] [g00] [g10] [g20] [g01] [g11] [g21] [g02] [g12] [g22] [m0] [m1] [m2] [c00] [c01] [c02] [c11] [c12] [c22] [d00] [d01] [d02] [d11] [d12] [d22] J_gamma
```
where the parameters follow the names given in *Eur.Phys.J.C 82 (2022) 1, 80*.
Please note, this version is designed to have parameters in both the numerator N and the denominator D floating in your fit. If you would like to fix parameters (e.g. in the denominator), then it can be (much) faster to calculate parts of the amplitude in the *calcUserVars()* function instead.

### Acknowledgement
This work is supported by the European Union’s Horizon 2020 research and innovation programme under Marie Sklodowska-Curie grant agreement No. 894790.
