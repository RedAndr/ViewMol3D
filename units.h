// Misc units
const double Clight=2.99792458e8  ;   // speed of light in m/s
const double Kboltz=3.166830e-6   ;   // Boltzmann constant
const double e2 = 14.399          ;   // Coulomb's law coeff if R in \AA and resulting E in eV
const double planck=6.6260755e-34 ;   // Planck's constant, in Js

// Distance units
const double bohr2ang = 0.5291772108 ;  // Conversion of length from bohr to angstrom
const double ang2bohr = 1/bohr2ang   ;

// Energy units
const double hartree2kcal = 627.5095; // Hartree to kcal/mol conversion
const double kcal2hartree = 1/hartree2kcal;

const double ev2kcal = 23.06054923 ;   // Conversion of energy in eV to energy in kcal/mol
const double kcal2ev = 1/ev2kcal   ;

const double Navog     = 6.022141510e23   ; // mol-1    Avogadro number
const double cal2joule = 4.184            ; // J/cal    calories to joule
const double ev2joule  = 1.6021765314e-19 ; // J/eV     electron volt to joule

const double hartree2joule = 4.35974417e-18  ; // Hatree to Joule conversion factor
const double joule2hartree = 1/hartree2joule;

// Mass units
const double amu2me = 1822.882   ;    // Conversion from mass in amu to mass in au (m_e)
const double me2amu = 1/amu2me   ;    // Conversion from mass in au (m_e) to mass in amu 


// Time units
const double tau2ps = 41341.447  ;    // Conversion from time in au to time in ps
const double ps2tau = 1/tau2ps   ;    // inverse

// Derived quantities
const double Rgas = Kboltz*hartree2kcal*1000.0; // gas constant R = 1.98722 cal/mole/K

const double pi=3.14159265358979323846;
