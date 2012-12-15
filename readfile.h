
#include <vector>

using namespace std;

struct atom { int type; double x,y,z,c; };  // atomic number, cartesian coordinates and electronic charge of the atom
struct bond { int a,b;  double o; };        // numbers of connected atoms and bond order = bond width
struct freq { int img;  double freqv,inten; vector<struct atom> coords; };
struct grad { int znuc; double x,y,z; };

typedef vector<struct atom> atomvector;
typedef vector<struct bond> bondvector;
typedef vector<struct freq> freqvector;
typedef vector<struct grad> gradvector;


class ReadMolecule {

  public:

    int
      tstrn,                  // total count of structures in file
      yes[4],                 // optimization flags in the Gaussian
      ircsave,                // save irc in file or not
      freqss,                 // is there frequences vectors or not
      gradss,                 // is there gradient vectors or not
      bread,                  // flag indicating is there bonds read or not
      cread;                  // coordinates read or not

    double 
      GraphData[512][2], 
      CurrentEnergy, 
      CurrentIRC;
    int NumData, CurrentPoint;

    int    ObjCount;
    double ObjData[256][32];

    atomvector atoms;
    bondvector bonds;
    freqvector freqs;
    gradvector grads;

    ReadMolecule();           // default constructor

    int OutFile (char* fname, int strn);
    int ALCHEMY (char* fname);
    int XYZ     (char* fname, int strn);
    void CenterMolecule();
    void ZOrientMolecule(int a1, int a2, int a3);
    void CalcPricipalAxes();

};

