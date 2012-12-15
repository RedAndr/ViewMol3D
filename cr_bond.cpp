
#include <math.h>
#include <vector>

using namespace std;

#include "newmatap.h"                   // newmat headers including advanced functions
#include "atom_rads.h"
#include "cr_bond.h"

#define len(x,y,z,xx,yy,zz) sqrt(sqr((x)-(xx))+sqr((y)-(yy))+sqr((z)-(zz)))
#define sqr(a)  ((a)*(a))

ReadMoleculeBond::ReadMoleculeBond(){
  BondRad0 =  0.08; 
  RadTol   = -1.40;
}

// To create/delete a bond
void ReadMoleculeBond::CreateBond(int a, int b) {

  // try to find bond to delete
  for ( vector <struct bond>::iterator it = bonds.begin(); 
        it!=bonds.end(); it++ ) {
     if ((( (*it).a==a+1) && ( (*it).b==b+1)) || 
         (( (*it).a==b+1) && ( (*it).b==a+1))) {
        bonds.erase( it );                                          // delete bond 
        return;                                                     // and exit    
     }                                                              
  }

  // if it was not found create it
  struct bond bond;

  if ( a > b ) { bond.a = a+1; bond.b = b+1; }
  else         { bond.a = b+1; bond.b = a+1; };
  double l = len(atoms[a].x,atoms[a].y,atoms[a].z,   atoms[b].x,atoms[b].y,atoms[b].z);
  double t = at_radv[atoms[a].type-1] + at_radv[atoms[b].type-1] + RadTol;
  bond.o = BondRad0 * t / l;
  bonds.push_back( bond );                                          // add new bond

}

#define NumAtoms (atoms.size())

// To create all bonds using density matrix from EHM
void ReadMoleculeBond::CreateBond(mindo3 *mcalc, atomvector_ehm atoms_ehm){

  int    i,j;
  struct bond bond;
  int   maxtype=0;

  for ( i=0; i<atoms.size(); i++ ) {
    if ( atoms[i].type > maxtype ) maxtype = atoms[i].type;
  }

  vector <struct atom>::iterator atom1, atom2;

  if ( (NumAtoms > 100) || (maxtype>10) ) {                                     // using van der Waals radii
     for ( i=1, atom1 = atoms.begin(); atom1!=atoms.end(); atom1++, i++ ) {
         for ( j=i+1, atom2 = atom1+1; atom2!=atoms.end(); atom2++, j++ ) {
            double l = len((*atom1).x,(*atom1).y,(*atom1).z,   (*atom2).x,(*atom2).y,(*atom2).z);
            double t = at_radv[(*atom1).type-1] + at_radv[(*atom2).type-1] + RadTol;
            if( l < t ) {
               if ( i > j ) { bond.a = i; bond.b = j; }
               else         { bond.a = j; bond.b = i; };
               bond.o = BondRad0 * t / l;
               bonds.push_back( bond );                                         // add new bond
            }
         }
     }

  } else {                                                                      // using Extended Huckel Method

     struct atom_ehm at;

     atoms_ehm.erase( atoms_ehm.begin(), atoms_ehm.end() );
     for (vector <struct atom>::iterator atit = atoms.begin(); atit!=atoms.end(); atit++) {
       at.atno =   (*atit).type;
       at.xyz[0] = (*atit).x;
       at.xyz[1] = (*atit).y;
       at.xyz[2] = (*atit).z;
       atoms_ehm.push_back(at);
     }
     
     mcalc = new mindo3(&atoms_ehm);
     mcalc->ehm::Energy();
     mcalc->calc_bondorders();

     for ( i=0; i<NumAtoms; i++ ) {
        for ( j=i+1; j<NumAtoms; j++ ) {
           double o = mcalc->get_bondorders( i, j );
           if ( o > 0.3 ) {
              bond.a = i+1; bond.b = j+1;
              bond.o = o/10;
              bonds.push_back( bond );                                          // add new bond
           }
        }
     }

//     delete mcalc;

  }

}
