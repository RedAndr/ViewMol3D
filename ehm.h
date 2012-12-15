
#include <vector>

#ifdef WIN32
#include <windows.h>
#endif

#include "newmatap.h"                   // newmat headers including advanced functions
#include "units.h"

#define  BOOL int
#define  TRUE 1
#define  FALSE 0

using namespace std;

// basis function structure
struct bfn {
  int    natom,                         // atom number
         type;                          // AO type: s,px,py or pz
  double ip;							// ionization potential
  double expn[6],                       // basis set
         coef[6];
};

// atom structure
struct atom_ehm { 
  int atno,Z;                           // atomic number and charge
  double xyz[3];                        // cartesian coordinate
  int nbf;                              // number of basis functions
};

typedef vector <struct atom_ehm> atomvector_ehm;

class ehm {

  public:

    atomvector_ehm *atoms;
    vector <bfn> bfns;

    int    nbf, nat,
           nclosed, MO_count;           // count of molecular orbitals
    BOOL   D_initialized;
    
    DiagonalMatrix Orbe; 
    SquareMatrix Orbs;
    SymmetricMatrix D,Orders;

    ehm (atomvector_ehm *_atoms);
    ~ehm ();

    int    numel (int charge);          // count of electrons included in calculation
    double Energy ();                   // calculated EHM energy

    void   print_mat(Matrix A);
    void   print_orbs ();
    void   print_dens ();

    double molorb(int numorb, double xyz[3]);
    double MOEnergy(int imo);
    double ElectronDensity(double xyz[3]);

    void   calc_bondorders();
    void   print_bondorders();
    double get_bondorders(int, int);

    double overlap      (bfn bfi, bfn bfj); 
    double doverlap     (bfn bfi, bfn bfj, int dir);
    double doverlap_num2(bfn bfi, bfn bfj, int dir);
    double doverlap_num4(bfn bfi, bfn bfj, int dir);

    double power(double x, double y);
    double dist2(double A[3],double B[3]);

    void   mkdens  ();

    double gto(bfn bf, double xyz[3]);

  private:

    void   init_basisfunctions ();

};


// orbitals names
const char or_nam[4][3] = {"S ","PX","PY","PZ"};

//STO-3G
//Robert F. Stewart, Small Gaussian Expansions of Slater-Type Orbitals, JCP 52 (1970), 431-438
//http://shay.ecn.purdue.edu/~barmstro/Research/pract/Chemistry/list/split/bassto_UM.ser.list

//#define stong 3
#define stong 6

#if stong==3

#define None      {                  0 ,                   0 ,                  0   }

#define gexps_1s  {    2.227660584e+00 ,     4.057711562e-01 ,     1.098175104e-01  }
#define gcoefs_1s {    1.543289673e-01 ,     5.353281423e-01 ,     4.446345422e-01  }

#define gexps_2s  {    2.581578398e+00 ,     1.567622104e-01 ,     6.018332272e-02  }
#define gcoefs_2s {   -5.994474934e-02 ,     5.960385398e-01 ,     4.581786291e-01  }

#define gexps_2p  {    9.192379002e-01 ,     2.359194503e-01 ,     8.009805746e-02  }
#define gcoefs_2p {    1.623948553e-01 ,     5.661708862e-01 ,     4.223071752e-01  }
                   
#define gexps_3s  {    5.641487709e-01 ,     6.924421391e-02 ,     3.269529097e-02  }
#define gcoefs_3s {   -1.782577972e-01 ,     8.612761663e-01 ,     2.261841969e-01  }
                   
#define gexps_3p  {    2.692880368e+00 ,     1.489359592e-01 ,     5.739585040e-02  }
#define gcoefs_3p {   -1.061945788e-02 ,     5.218564264e-01 ,     5.450015143e-01  }

#else

//STO-6G
//Robert F. Stewart, Small Gaussian Expansions of Slater-Type Orbitals, JCP 52 (1970), 431-438

#define None      {0,                0,                0,               0,              0,              0}

#define gexps_1s  {2.310303149e01   ,4.235915534e00   ,1.185056519e00  ,4.070988982e-01,1.580884151e-01,6.510953954e-02}
#define gcoefs_1s {9.163596280e-03  ,4.936149294e-02  ,1.685383049e-01 ,3.705627997e-01,4.164915298e-01,1.303340841e-01}

#define gexps_2s  {2.768496241e01   ,5.077140627e00   ,1.426786050e00  ,2.040335729e-01,9.260298399e-02,4.416183978e-02}
#define gcoefs_2s {-4.151277819e-03 ,-2.067024148e-02 ,-5.150303337e-02,3.346271174e-01,5.621061301e-01,1.712994697e-01}

#define gexps_2p  {5.868285913e00   ,1.530329631e00   ,5.475665231e-01 ,2.288932733e-01,1.046655969e-01,4.948220127e-02}
#define gcoefs_2p {7.924233646e-03  ,5.144104825e-02  ,1.898400060e-01 ,4.049863191e-01,4.012362861e-01,1.051855189e-01}

#define gexps_3s  {3.273031938e00   ,9.200611311e-01  ,3.593349765e-01 ,8.636686991e-02,4.797373812e-02,2.724741144e-02}
#define gcoefs_3s {-6.775596947e-03 ,-5.639325779e-02 ,-1.587856086e-01,5.534527651e-01,5.015351020e-01,7.223633674e-02}

#define gexps_3p  {5.077973607e00   ,1.340786940e00   ,2.248434849e-01 ,1.131741848e-01,6.076408893e-02,3.315424265e-02}
#define gcoefs_3p {-3.329929840e-03 ,-1.419488340e-02 ,1.639395770e-01 ,4.485358256e-01,3.908813050e-01,7.411456232e-02}

#endif


const int NQN[19] = {  0, 1, 1,                         // principle quantum number N
        2, 2, 2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3  };


const double gexps[3][2][stong] = {                     // indexed by N,s_or_p:
  { gexps_1s, None     },  // N=1
  { gexps_2s, gexps_2p },  // N=2
  { gexps_3s, gexps_3p },  // N=3
}; //  s       p
    
    
const double gcoefs[3][2][stong] = {                    // indexed by N,s_or_p:
  { gcoefs_1s, None      },  // N=1
  { gcoefs_2s, gcoefs_2p },  // N=2
  { gcoefs_3s, gcoefs_3p },  // N=3
}; //  s       p

const int s_or_p[4] = {0,1,1,1};                        // whether the func is s or p type, based on the L QN


// s atomic orbital ionization potential for two center resonance integral term
const double IPs[19] = { 0.0, -13.605, -24.6,
        -5.4, -10.0, -15.160, -21.340, -27.510, -35.300, -43.700, -17.820,
         0.0,   0.0,   0.0  , -17.82,  -21.100, -23.840, -25.260, 0.0};

// p atomic orbital ionization potential for two center resonance integral term
const double IPp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, -8.520, -11.540, -14.340, -17.910, -20.890, -8.510,
        0.0, 0.0,  0.0,    -8.51,  -10.290, -12.410, -15.090,  0.0};

// s-type Slater atomic orbital exponent
const double zetas[19] = { 0.0, 
          1.30,                                                 2.0925,
          0.650, 0.975, 1.211156, 1.739391, 2.704546, 3.640575, 3.111270, 0.0,
          0.0,   0.0,   0.0,      1.629173, 1.926108, 1.719480, 3.430887, 0.0};

// p-type Slater atomic orbital exponent
const double zetap[19] = { 0.0, 
          0.0,                                                        0.0,
          0.0, 0.0, 0.972826, 1.709645, 1.870839, 2.168448, 1.419860, 0.0,
          0.0, 0.0, 0.0,      1.381721, 1.590665, 1.403205, 1.627017, 0.0};

const int nbfat[19] = { 0, 
          1,                   1,
          0, 0, 4, 4, 4, 4, 4, 0,
          0, 0, 0, 4, 4, 4, 4, 0};

const int CoreQ[19] = { 0, 
          1,                   1,
          0, 0, 3, 4, 5, 6, 7, 0,
          0, 0, 0, 4, 5, 6, 7, 0};
