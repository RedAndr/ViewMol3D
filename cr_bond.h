
#include "readfile.h"
#include "mindo3.h"

class ReadMoleculeBond : public ReadMolecule {

  double BondRad0; 
  double RadTol;
  
  public:
  
  ReadMoleculeBond();
  void CreateBond(int a, int b);
  void CreateBond(mindo3 *mcalc, atomvector_ehm atoms_ehm);

};

