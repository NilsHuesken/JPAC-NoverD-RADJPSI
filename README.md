# JPAC-NoverD-RADJPSI

This code reproduces the results shown in *A. Rodas et al. (JPAC Collaboration), Eur.Phys.J.C 82 (2022) 1, 80 (DOI: 10.1140/epjc/s10052-022-10014-8, arXiv:2110.00027*), using the parameter values given in the Appendix.
The dispersion integrals in eq. (4) for the different choices of rhoN have been provided in the form of lookup-tables by Alessandro Pilloni and Arkaitz Rodas.

For matrix calculations, the code uses the *Eigen library* (see https://eigen.tuxfamily.org) - a copy of Eigen is **not** included in this repository. To run the code, you need to download a copy of Eigen and adjust the include-path for Eigen/Dense.

ROOT is used to plot the squared-magnitude of the amplitude compared to binned BESIII data for the pipi- (*Phys.Rev.D 92 (2015) 5, 052003*) and KK-channels (*Phys.Rev.D 98 (2018) 7, 072003*).

If useful, please feel free to try and use this amplitude as a description of the pipi or KK S- and D-wave between 1 and 2.5 GeV in your work, citing the original publication *Eur.Phys.J.C 82 (2022) 1, 80* as well as https://eigen.tuxfamily.org.
You can also refer to doi.org/XXXX/zenodo.XXXX for this repository.

### Acknowledgement
This work is supported by the European Union’s Horizon 2020 research and innovation programme under Marie Sklodowska-Curie grant agreement No. 894790.
