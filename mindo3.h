
#include <vector>
#include "ehm.h"

using namespace std;

class mindo3 : public ehm {

  public:

    int    MaxSCF, SCFit;
    double TolSCF, Ediff;
    BOOL   F0_initialized;
    
    mindo3 (atomvector_ehm *_atoms);
    ~mindo3();

    int    numel (int charge);
    double Energy ();

    void   forces ( Matrix *forces );
    void   num_forces_right   ( Matrix *forces, double dx );
    void   num_forces_left    ( Matrix *forces, double dx );
    void   num_forces_central ( Matrix *forces, double dx );
    void   num_forces_central4( Matrix *forces, double dx );
    void   num_hessian_right  ( SymmetricMatrix *hess, double dx );
    void   num_hessian_central( SymmetricMatrix *hess, double dx );

    void   print_grads( Matrix G);
    void   print_hess ( SymmetricMatrix H);

    void   ElectronDensity_renormalize();

  private:

    SymmetricMatrix F0,F1,F2,F;

    double axy[19][19], Bxy[19][19];

    void   init_parameters ();

    double gamma(struct atom_ehm ati, struct atom_ehm atj);
    double scale(int atnoi,int atnoj,double R);

    double enuke ();
    double refeng();

    double g ( bfn bfi, bfn bfj );
    double h ( bfn bfi, bfn bfj );

    void   calc_F0 ();
    void   calc_F1 ();
    void   calc_F2 ();

    void   guess_D ();
    double SCF     ();

};

// MINDO/3 Parameters: Thru Ar in eV

//s and p atomic orbital one-electron one-center integrals
const double Uss[19] = { 0.0, -12.505, 0.0,
        0.0, 0.0, -33.61, -51.79, -66.06, -91.73, -129.86, 0.0,
        0.0, 0.0,   0.0 , -39.82, -56.23, -73.39,  -98.99, 0.0}; 
const double Upp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, -25.11, -39.18, -56.40, -78.80, -105.93, 0.0,
        0.0, 0.0,   0.0 , -29.15, -42.31, -57.25,  -76.43, 0.0};

// s-s atomic orbital one center two electron repulsion integral
const double gss[19] = { 0.0, 12.848, 0.0,
        0.0, 0.0, 10.59, 12.23, 13.59, 15.42, 16.92, 0.0,
        0.0, 0.0,  0.0 ,  9.82, 11.56, 12.88, 15.03, 0.0};

// s-p atomic orbital one center two electron repulsion integral
const double gsp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, 9.56, 11.47, 12.66, 14.48, 17.25, 0.0,
        0.0, 0.0,  0.0,  8.36, 10.08, 11.26, 13.16, 0.0};

// p-p atomic orbital one center two electron repulsion integral
const double gpp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, 8.86, 11.08, 12.98, 14.52, 16.71, 0.0,
        0.0, 0.0, 0.0 ,  7.31,  8.64,  9.90, 11.30, 0.0};

// p-p' atomic orbital one center two electron repulsion integral
const double gppp[19] = { 0.0, 0.0, 0.0,
         0.0, 0.0, 7.86, 9.84, 11.59, 12.98, 14.91, 0.0,
         0.0, 0.0,  0.0, 6.54,  7.68,  8.83,  9.97, 0.0};

// s-p atomic orbital one-center two-electron exchange integral
const double hsp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, 1.81, 2.43, 3.14, 3.94, 4.83, 0.0,
        0.0, 0.0, 0.0,  1.32, 1.92, 2.26, 2.42, 0.0};

const double hppp[19] = { 0.0, 0.0, 0.0,
         0.0, 0.0, 0.50, 0.62, 0.70, 0.77, 0.90, 0.0,
         0.0, 0.0, 0.0 , 0.38, 0.48, 0.54, 0.67, 0.0};

// averaged repulsion integral for use in gamma
const double f03[19] = { 0.0, 
        12.848,                                            10.0,                             
        10.0,  0.0, 8.958, 10.833, 12.377, 13.985, 16.250, 10.0, 
        10.0,  0.0, 0.000,   7.57,  9.00 , 10.20 , 11.73,  10.0};

// Atomic heat of formations: Mopac got from CRC
const double Hfat[19] = { 0.0, 
         52.102, 0.0,
         0.0, 0.0, 135.7, 170.89, 113.0, 59.559, 18.86, 0.0,
         0.0, 0.0, 0.0, 106.0, 79.8, 65.65, 28.95, 0.0};

// Default isolated atomic energy values from Mopac:EISOL3
const double Eat[19] = {0.0, 
       -12.505, 0.0,
       0.0 ,0.0,-61.70,-119.47,-187.51,-307.07,-475.00,0.0,
       0.0,0.0,0.0,-90.98,-150.81,-229.15,-345.93,0.0};

