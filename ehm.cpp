//
// Extended Huckel method
// STO-3G or STO-6G orbitals
//

#include <cmath>
#include "atom_rads.h"
#include "ehm.h"
#include "debug.h"

double ehm::power(double x, double y) {
  if ((x==0) && (y==0)) return 1.0;
  if  (x==0)            return 0.0;
  if  (x==1)            return 1.0;
  if  (y==2)            return x*x;
  if  (y==3)            return x*x*x;
  if  (y==0.5)          return sqrt(x);
  if  (y==-0.5)         return 1/sqrt(x);
  return pow(fabs(x),y);
}


ehm::ehm(atomvector_ehm *_atoms) : atoms(_atoms) {

  nat = atoms->size();

  init_basisfunctions();

  D_initialized  = FALSE;

  nclosed  = numel(0)/2;                                                       // number of occupied orbitals, closed shells only
  MO_count = bfns.size();                                                      // number of MO
  nbf      = bfns.size();                                                      // number of basis functions in an atom_ehm list;

}

ehm::~ehm() {

  Orbe.release(); 
  Orbs.release();
  D   .release();
  Orders.release();

}


void ehm::init_basisfunctions() {

  bfns.erase( bfns.begin(), bfns.end() );

  for ( int iat=0, na=0; iat<nat; iat++, na++) {

     atom_ehm *at = &((*atoms)[iat]);
     int atno = (*at).atno;
          
     (*at).Z    = CoreQ [ atno ];
     (*at).nbf  = nbfat [ atno ];

     for (int i=0; i<(*at).nbf; i++) {
        struct bfn bfn;
        double zeta;
        bfn.natom = na;
        bfn.type  = i;
        if ( i==0 ) {
            zeta   = zetas[atno];
            bfn.ip = IPs  [atno];
        } else {
            zeta   = zetap[atno];
            bfn.ip = IPp  [atno];
        }
        for (int j=0;j<stong;j++) {
            bfn.expn[j] = gexps  [ NQN[atno]-1 ] [ s_or_p[i] ] [ j ] * zeta * zeta;
            bfn.coef[j] = gcoefs [ NQN[atno]-1 ] [ s_or_p[i] ] [ j ];
        }
        bfns.push_back(bfn);
     }

  }

}



double ehm::dist2(double A[3],double B[3]) {
    return power(A[0]-B[0],2)+power(A[1]-B[1],2)+power(A[2]-B[2],2);
}


// Number of electrons in the system
int ehm::numel(int charge) {

    int nel = 0;

    for ( int i=0; i<nat; i++ ) {
      nel += (*atoms)[i].Z;
    }

    return nel-charge;

}


// Calculate the overlap integrals using a gaussian expansion STO-nG
double ehm::overlap(bfn bfi, bfn bfj) { 

    int i,j;
    double ri[3], rj[3];                                                        // distance in bohr

    for ( i=0; i<3; i++ ) {
      ri[i]=(*atoms)[bfi.natom].xyz[i]/bohr2ang;
      rj[i]=(*atoms)[bfj.natom].xyz[i]/bohr2ang;
    }

    double RR = power(ri[0]-rj[0],2) + power(ri[1]-rj[1],2) + power(ri[2]-rj[2],2);
    int itype = bfi.type;
    int jtype = bfj.type;
    double Sij = 0.0;

    for ( i=0; i<stong; i++) {
      for ( j=0; j<stong; j++) {

            double amb = bfi.expn[i] + bfj.expn[j];
            double apb = bfi.expn[i] * bfj.expn[j];
            double adb = apb/amb;                             
            double tomb, abn;

            if ((itype > 0) && (jtype > 0)) {                                   // is = 4
                tomb = (ri[itype-1]-rj[itype-1])*(ri[jtype-1]-rj[jtype-1]);
                abn = -adb*tomb;
                if (itype == jtype) abn += 0.5;
                abn = 4*abn*sqrt(apb)/amb;
            } else 
            
            if (itype > 0) {                                                    // is = 3
                tomb = (ri[itype-1]-rj[itype-1]);
                abn = -2*tomb*bfj.expn[j] * sqrt(bfi.expn[i])/amb;
            } else 
            
            if (jtype > 0) {                                                    // is = 2
                tomb = (ri[jtype-1]-rj[jtype-1]);
                abn =  2*tomb*bfi.expn[i] * sqrt(bfj.expn[j])/amb;
            } else {
 
                abn = 1.0;                                                      // is = 1
            }    
            
            if (adb*RR < 90) {
                Sij += bfi.coef[i] * bfj.coef[j] * pow(2*sqrt(apb)/amb,1.5) * exp(-adb*RR)*abn;
            }
      }
    }

    return Sij;

}


double ehm::doverlap_num2(struct bfn bfi, struct bfn bfj, int dir) {

   double eps = 1e-6;

   (*atoms)[bfi.natom].xyz[dir] += eps;  
   double ovr1 = overlap(bfi,bfj);
   (*atoms)[bfi.natom].xyz[dir] -= eps;  

   (*atoms)[bfi.natom].xyz[dir] -= eps;  
   double ovr2 = overlap(bfi,bfj);
   (*atoms)[bfi.natom].xyz[dir] += eps;

   return (ovr1-ovr2)/(eps+eps);

}


double ehm::doverlap_num4(struct bfn bfi, struct bfn bfj, int dir) {

   double eps = 1e-6;

   (*atoms)[bfj.natom].xyz[dir] += 2*eps;  
   double ovr1 = overlap(bfi,bfj);
   (*atoms)[bfj.natom].xyz[dir] -= 2*eps;  

   (*atoms)[bfj.natom].xyz[dir] += eps;  
   double ovr2 = overlap(bfi,bfj);
   (*atoms)[bfj.natom].xyz[dir] -= eps;  

   (*atoms)[bfj.natom].xyz[dir] -= eps;  
   double ovr3 = overlap(bfi,bfj);
   (*atoms)[bfj.natom].xyz[dir] += eps;

   (*atoms)[bfj.natom].xyz[dir] -= 2*eps;  
   double ovr4 = overlap(bfi,bfj);
   (*atoms)[bfj.natom].xyz[dir] += 2*eps;

   return (ovr1-8*ovr2+8*ovr3-ovr4)/(12*eps);

}

// Analytical doverlap function
double ehm::doverlap(struct bfn bfi, struct bfn bfj, int dir) {

    int i,j;
    double DS = 0.0, SS, amb, apb, adb, abn, adr, del1, del2 ,del3;

    double ri[3], rj[3];                                                        // distance in bohr
    for ( i=0; i<3; i++ ) {
      ri[i]=(*atoms)[bfi.natom].xyz[i] * ang2bohr;
      rj[i]=(*atoms)[bfj.natom].xyz[i] * ang2bohr;
    }
    double RR = power(ri[0]-rj[0],2) + power(ri[1]-rj[1],2) + power(ri[2]-rj[2],2);
    del1 = ri[dir] - rj[dir];

    int itype = bfi.type;
    int jtype = bfj.type;

    for ( i=0; i<stong; i++ ) {                                 

      for ( j=0; j<stong; j++ ) {                               

            amb = bfi.expn[i] + bfj.expn[j];                               
            apb = bfi.expn[i] * bfj.expn[j]; 
            adb = apb/amb;
            adr = min(adb*RR,35.0);

            if ((itype == 0) && (jtype == 0)) {                                 // is = 1 (s|s)
                    abn = -2*adb*del1;
            } else if ((itype == 0) && (jtype > 0)) {                           // sp

                if ((jtype-1) == dir) {                                         // is = 3 (s/p)
                    abn = (2*adb/sqrt(bfj.expn[j]))*(1-2*adb*del1*del1);
                } else {                                                        // is = 2 (s/p')
                    del2 = ri[jtype-1]-rj[jtype-1];
                    abn = -4*adb*adb*del1*del2/sqrt(bfj.expn[j]);
                }
            } else if ((itype > 0) && (jtype == 0)) {                           // ps
                if ((itype-1) == dir) {                                         // is = 5 (p/s)
                    abn = -2*adb/sqrt(bfi.expn[i])*(1-2*adb*del1*del1);
                } else {                                                        // is = 4 (p'/s)
                    del2 = ri[itype-1]-rj[itype-1];
                    abn = 4*adb*adb*del1*del2/sqrt(bfi.expn[i]);
                }

            } else if (itype == jtype) {                                        // pp
                if (dir == (itype-1)) {                                         // is = 9 (p|p)
                    abn=-8*adb*adb*del1/sqrt(apb)*(1.5-adb*del1*del1);
                } else {                                                        // is = 8 (p'|p')
                    del2 = ri[jtype-1]-rj[jtype-1];
                    abn=-8*pow(adb,2)*del1/sqrt(apb)*(0.5-adb*del2*del2);
                }

            } else if ((dir != (itype-1)) && (dir != (jtype-1))) {              // is = 7(p'|p")
                del2 = ri[itype-1] - rj[itype-1];
                del3 = ri[jtype-1] - rj[jtype-1];
                abn=8*pow(adb,3)*del1*del2*del3/sqrt(apb);

            } else {                                                            // is = 6 (p|p') or (p'|p)
                del2 = ri[itype+jtype-dir-2]-rj[itype+jtype-dir-2];
                abn=-4*adb*adb*del2/sqrt(apb)*(1-2*adb*del1*del1);
            }

            SS = sqrt(power(2*sqrt(apb)/amb,3)) * exp(-adr) * abn;

            DS += SS * bfi.coef[i] * bfj.coef[j];

      }

    }

    return DS * ang2bohr;

}

void ehm::print_orbs() {

  char sss[256]; 

  sprintf(sss,"%4s %4s "," ","E=");DEBUG(sss);
  for ( int j=0; j<nbf; j++ ) {
    sprintf(sss, "%11.6lf",Orbe(j+1)/*/ev2kcal*/); 
    DEBUG(sss);
  }
  DEBUG("\n");
  for ( int i=0; i<nbf; i++ ) {
      sprintf(sss,"%4s %4s ",at_nam[(*atoms)[bfns[i].natom].atno-1], or_nam[bfns[i].type]);
      DEBUG(sss);
      for(int j=0;j<nbf;j++) {
        sprintf(sss,"%11.6lf",Orbs(i+1,j+1));
        DEBUG(sss);
      }
      DEBUG("\n");
  }
}


void ehm::print_dens() {

  char sss[256]; int i;

  sprintf(sss,"%4s %4s "," "," ");DEBUG(sss);
  for ( i=0; i<nbf; i++ ) {
      sprintf(sss," %4s %4s ",at_nam[(*atoms)[bfns[i].natom].atno-1], or_nam[bfns[i].type]);
      DEBUG(sss);
  }
  DEBUG("\n");
  for ( i=0; i<nbf; i++ ) {
      sprintf(sss,"%4s %4s ",at_nam[(*atoms)[bfns[i].natom].atno-1], or_nam[bfns[i].type]);
      DEBUG(sss);
      for(int j=0;j<nbf;j++) {
        if ( j > i ) {
          DEBUG("           ");
        } else {
          sprintf(sss,"%11.6lf",D(i+1,j+1));
          DEBUG(sss);
        }
      }
      DEBUG("\n");
  }
}

void ehm::print_mat(Matrix A) {

  char sss[256]; 

  for ( int i=0; i<nbf; i++ ) {
      for(int j=0;j<nbf; j++) {
        sprintf(sss,"%11.6lf",A(i+1,j+1));
        DEBUG(sss);
      }
      DEBUG("\n");
  }
}


// Gaussian type orbital, AO = GTO
double ehm::gto(struct bfn bf, double xyz[3]) {

  double R2 = dist2( (*atoms)[bf.natom].xyz, xyz );

  if ( R2 > 20.0) return 0.0; else {

    double sum = 0.0; 

    for ( int i=0; i<stong; i++ ) {
      sum += bf.coef[i] * exp( -bf.expn[i] * R2 );
    }

    if (bf.type == 0)
       return sum;                        // S
    else                
       return sum * xyz[ bf.type-1 ];     // Pxyz

  }

}


// Molecular orbital, MO = LCAO = LCGTO
double ehm::molorb(int numorb, double xyz[3]) {

  double sum = 0.0;

  for ( int i=0; i<nbf; i++ ) {
      sum += Orbs(i+1,numorb+1) * gto(bfns[i],xyz);
  }

  return sum;
}


// Energy of molecular orbital 
double ehm::MOEnergy(int imo) {
  return Orbe(imo+1);
}


// Calculate electron density for point
double ehm::ElectronDensity(double xyz[3]) {

  int i;  

  vector <double> gtos(nbf);
  for ( i=0; i<nbf; i++) {
    gtos[i] = gto(bfns[i],xyz);
  }

  double sum = 0;
  for ( i=0; i<nbf; i++) {
     for ( int j=0; j<nbf; j++ ) {
         sum += gtos[i] * gtos[j] * D(i+1, j+1);
     }
  }

  gtos.clear();

  return sum;

}


// Bond orders functions
void ehm::calc_bondorders(){

  Orders.resize(nat);
  Orders = 0;

  for ( int i=1; i<=nbf; i++ )
    for ( int j=i; j<=nbf; j++ )
       Orders ( bfns[i-1].natom+1, bfns[j-1].natom+1 ) += power( D(i,j), 2 );

}

// Print matrix of bond orders
void ehm::print_bondorders() {

  int i;
  
  printf(" %11s "," ");
  for ( i=0; i<nat; i++ ) {
    printf(" %8s%3u ",at_nam[(*atoms)[i].atno-1],i+1);
  }
  puts("");
  for ( i=1; i<=nat; i++ ) {
    printf(" %8s%3u ",at_nam[(*atoms)[i-1].atno-1],i);
    for ( int j=1; j<=nat; j++ )
      if ( i<j )
        printf(" %11s "," ");
      else
        printf("%13.4lf",Orders(i,j));
    puts("");
  }
}

// Get order for particular bond
double ehm::get_bondorders(int i, int j) {
  return Orders(i+1,j+1);
}

// Form a density matrix C*Ct given eigenvectors C[nstart:nstop,:]
void ehm::mkdens() {
    D << Orbs.columns(1,nclosed) * Orbs.columns(1,nclosed).t() * 2;
}

// Extended Huckel Method Energy
double ehm::Energy() {

  SymmetricMatrix S(nbf), H(nbf);
  DiagonalMatrix d(nbf);
  SquareMatrix P(nbf);

  double K = 1.75;
  int i,j;

  for ( i=1; i<=nbf; i++)
    for ( j=1; j<=nbf; j++)
      if (i>=j) {
        S(i,j) = overlap ( bfns[i-1], bfns[j-1] );                             // overlap matrix
        if (i!=j) {
          double I = bfns[i-1].ip + bfns[j-1].ip;                              // Ionization potential for an atoms couple
          H(i,j) = 0.5 * K * S(i,j) * I;                                       // Fock matrix
        } else {
          H(i,i) = bfns[i-1].ip;                                               // Ionization potential for an atom
        }
      }

  EigenValues( S, d, P );                                                       // Diagonalize S = P*D*P^(-1)
  for ( i=1; i<=nbf; i++) d(i,i) = power( d(i,i), -0.5);                        // D^(-1/2)
  S << P*d*P.t();                                                               // Ortogonalize S matrix => S^(-1/2)

  H << S*H*S;                                                                   // Transformation H to H'

  Orbe.resize(nbf);
  Orbs.resize(nbf);

  EigenValues( H, Orbe, Orbs );                                                 // Diagonalize H'
//  Orbs << S * Orbs;                                                           // Transformation C' to C

  S.release();
  H.release();
  d.release();
  P.release();

  mkdens();                                                                     // Calculate D
  D_initialized  = TRUE;

  double Eel=0;
  for ( i=1; i<=nclosed; i++) Eel += Orbe(i);
  Eel *= 2;
    
  return Eel;

}

