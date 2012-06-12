//===========================================================================================================
// Author: Jens Chluba
// first implementation: April 2010
// last modification: June 2012
// purpose: to access the main CosmoRec code
//===========================================================================================================
// 08.06.2012: Added version to communicate Hubble 

#ifndef COSMOREC_H
#define COSMOREC_H

//===========================================================================================================
// CosmoRec for run from console
//===========================================================================================================
int CosmoRec(int narg, char *args[]);

//===========================================================================================================
// to call CosmoRec with a list of parameters in batch mode (e.g. when calling it from CosmoMC).
// 
// runmode == 0: CosmoRec run with diffusion
//            1: CosmoRec run without diffusion
//            2: Recfast++ run (equivalent of the original Recfast version)
//            3: Recfast++ run with correction function of Chluba & Thomas, 2010
//
// On entry, the array z_arr should contain the redshifts at which Xe and Te are required. nz 
// determines the number of redshift points. Xe_arr and Te_arr will contain the solution on exit.
//
// Furthermore, runpars[0] defines the dark matter annihilation efficiency in eV/s.
// runpars[1] switches the accuracy of the recombination model:
//
// runpars[1]==-1: closest equivalent of 'HyRec' case (Haimoud & Hirata, 2010)
// runpars[1]== 0: default
// runpars[1]== 1: 2g for n<=4 & Raman for n<=3
// runpars[1]== 2: 2g for n<=8 & Raman for n<=7
// runpars[1]== 3: 2g for n<=8 & Raman for n<=7 + Helium feedback up to n=5
//
// The value of runpars[1] is only important for runmode 0 & 1.
//===========================================================================================================
int CosmoRec(const int runmode, const double runpars[5], 
             const double omegac, const double omegab, 
             const double omegak, const double Nnu,  
             const double h0, const double tcmb, const double yhe, 
             const int nz, double *z_arr, double *Xe_arr, double *Te_arr,
             const int label);


//===========================================================================================================
// Wrap the C++ Fortran routine to be allow calling from Fortran. Arguments are as above.
// Added 06.03.2011 (Richard Shaw)
//===========================================================================================================
// 08.06.2012: Added version to communicate Hubble
//===========================================================================================================
extern "C" {
    
    void cosmorec_calc_cpp_(const int * runmode, const double * runpars, 
                            const double * omega_c, const double * omega_b, const double * omega_k, 
                            const double * num_nu, const double * h0, 
                            const double * t_cmb, const double * y_he, 
                            double * za_in, double * xe_out, double * tb_out, 
                            const int * len, const int* label);
    
    void cosmorec_calc_h_cpp_(const int * runmode, const double * runpars, 
                              const double * omega_c, const double * omega_b, const double * omega_k, 
                              const double * num_nu, const double * h0, 
                              const double * t_cmb, const double * y_he, 
                              const double * z_Hz, const double * Hz, const int * nz,
                              double * za_in, double * xe_out, double * tb_out, 
                              const int * len, const int* label); 
}

//===========================================================================================================
// to cleanup after finishing with CosmoRec
//===========================================================================================================
void cleanup_CosmoRec();

//===========================================================================================================
#endif