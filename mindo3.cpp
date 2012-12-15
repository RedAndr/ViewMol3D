// Modified Intermediate Neglect of Differential Overlap, version 3 (MINDO/3)
// semiempirical method by Dewar, reference:
// Bingham, R. C., Dewar, M. J. S. and Lo, D. H. JACS, 97, 1285, 1307, (1975)
// See also information about SM methods:
// http://www.cachesoftware.com/mopac/Mopac2002manual/node439.html

#include <cmath>
#include "mindo3.h"
//#include "debug.h"


mindo3::mindo3(atomvector_ehm *_atoms) : ehm(_atoms) {

  init_parameters();

  MaxSCF = 256/*64*/;
  TolSCF = 1e-9;//1.0e-7;

  F0_initialized = FALSE;

}


mindo3::~mindo3 () {

  F0.release();
  F1.release();
  F2.release();
  F .release();

}


// Si parameters from Edwards, Arthur H.; Fowler, W. Beall.  
// Semiempirical molecular orbital techniques applied to silicon dioxide:  MINDO/3. 
// Journal of Physics and Chemistry of Solids  (1985),  46(7),  841-57
void mindo3::init_parameters() {
  int i, j;

  for(i=0;i<19;i++) for(j=0;j<19;j++) axy[i][j] = 2.0;		// fill by default

// axy Core repulsion function terms
  axy[1][1 ] = 1.489450; axy[1][5 ] = 2.090352; axy[1][6 ] = 1.475836; axy[1][7 ] = 0.589380; axy[1][8 ] = 0.478901; 
  axy[1][9 ] = 3.771362; axy[1][14] = 0.940789; axy[1][15] = 0.923170; axy[1][16] = 1.700689; axy[1][17] = 2.089404; 
  axy[5][5 ] = 2.280544; axy[5][6 ] = 2.138291; axy[5][7 ] = 1.909763; axy[5][8 ] = 2.484827; axy[5][9 ] = 2.862183; 
  axy[6][6 ] = 1.371208; axy[6][7 ] = 1.635259; axy[6][8 ] = 1.820975; axy[6][9 ] = 2.725913; axy[6][14] = 1.101382;  
  axy[6][15] = 1.029693; axy[6][16] = 1.761370; axy[6][17] = 1.676222; axy[7][7 ] = 2.209618; axy[7][8 ] = 1.873859; 
  axy[7][9 ] = 2.861667; axy[8][8 ] = 1.537190; axy[8][9 ] = 2.266949; axy[8][14] = 1.3716248;axy[9][9 ] = 3.864997; 
  axy[14][14]= 0.918432; axy[15][15]= 1.186652; axy[16][16]= 1.751617; axy[17][17]= 1.792125;

  axy[ 8][14] = 1.5985778157;
  axy[14][14] = 0.6620676282;   
         
  for(i=0;i<19;i++) for(j=i+1;j<19;j++) axy[j][i] = axy[i][j];


  for(i=0;i<19;i++) for(j=0;j<19;j++) Bxy[i][j] = 0.5;		// fill by default

// Diatomic two center one-electron resonance integral multiplier
  Bxy[1 ][1 ] = 0.244770;  Bxy[1 ][5 ] = 0.185347;  Bxy[1 ][6 ] = 0.315011;  Bxy[1 ][7 ] = 0.360776;  Bxy[1 ][8 ] = 0.417759; 
  Bxy[1 ][9 ] = 0.195242;  Bxy[1 ][14] = 0.289647;  Bxy[1 ][15] = 0.320118;  Bxy[1 ][16] = 0.220654;  Bxy[1 ][17] = 0.231653;
  Bxy[5 ][5 ] = 0.151324;  Bxy[5 ][6 ] = 0.250031;  Bxy[5 ][7 ] = 0.310959;  Bxy[5 ][8 ] = 0.349745;  Bxy[5 ][9 ] = 0.219591;
  Bxy[6 ][6 ] = 0.419907;  Bxy[6 ][7 ] = 0.410886;  Bxy[6 ][8 ] = 0.464514;  Bxy[6 ][9 ] = 0.247494;  Bxy[6 ][14] = 0.411377; 
  Bxy[6 ][15] = 0.457816;  Bxy[6 ][16] = 0.284620;  Bxy[6 ][17] = 0.315480;  Bxy[7 ][7 ] = 0.377342;  Bxy[7 ][8 ] = 0.458110; 
  Bxy[7 ][9 ] = 0.205347;  Bxy[8 ][8 ] = 0.659407;  Bxy[8 ][9 ] = 0.334044;  Bxy[8 ][14] = 0.5701762; Bxy[9 ][9 ] = 0.197464;  
  Bxy[14][14] = 0.291703;  Bxy[15][15] = 0.311790;  Bxy[16][16] = 0.202489;  Bxy[17][17] = 0.258969;

  Bxy[ 8][14] = 0.4595320949;   
  Bxy[14][14] = 0.3836309113;

  for(i=0;i<19;i++) for(j=i+1;j<19;j++) Bxy[j][i] = Bxy[i][j];

}      


// Coulomb repulsion that goes to the proper limit at R=0"
double mindo3::gamma(struct atom_ehm ati, struct atom_ehm atj) {
    return e2/sqrt( dist2( ati.xyz, atj.xyz ) + 0.25*power(e2/f03[ati.atno] + e2/f03[atj.atno],2) );
}


// Prefactor from the nuclear repulsion term
double mindo3::scale(int atnoi, int atnoj, double R) {

    double alpha = axy[atnoi][atnoj];                                   // Part of the scale factor for the nuclear repulsion

    if (atnoi == 1) {
        if ((atnoj == 7) || (atnoj == 8))
            return alpha*exp(-R);
    } else 
      if (atnoj == 1)
        if ((atnoi == 7) || (atnoi == 8))
            return alpha*exp(-R);

    return exp(-alpha*R);

}


// Compute the nuclear repulsion energy
double mindo3::enuke() {
    double enuke = 0.0;
    for ( int i=0; i<nat; i++ ) {
      struct atom_ehm ati = (*atoms)[i];
      for ( int j=i+1; j<nat; j++ ) {
            struct atom_ehm atj = (*atoms)[j];
            double R = sqrt( dist2( ati.xyz, atj.xyz ) );
            double sc = scale( ati.atno, atj.atno, R );
            double gammaij = gamma( ati, atj );
            enuke += ati.Z*atj.Z*(gammaij + (e2/R-gammaij)*sc);
       }
    }
    return enuke;
}


// Ref = heat of formation - energy of atomization
double mindo3::refeng() {

    double eat = 0.0,
          hfat = 0.0;

    for ( int i=0; i<nat; i++ ) {
        eat  += Eat [(*atoms)[i].atno];
        hfat += Hfat[(*atoms)[i].atno];
    }
    return hfat-eat*ev2kcal;
}


// Form the zero-iteration (density matrix independent) Fock matrix
void mindo3::calc_F0() {

  if ( F0_initialized ) return;

  F0 = 0;

  for ( int i=1; i<=nbf; i++ ) {

     int iat = bfns[i-1].natom;
//     F0(i,i) = bfns[i-1].u;
     if ( bfns[i-1].type==0 )  F0(i,i) = Uss [ (*atoms)[iat].atno ];
     else                      F0(i,i) = Upp [ (*atoms)[iat].atno ];
     int jato = -1;

     for ( int j=1; j<=nbf; j++ ) {

        int jat = bfns[j-1].natom;

        if ( iat != jat ) {

           if ( jat != jato ) { 
              F0(i,i) -= gamma( (*atoms)[iat], (*atoms)[jat] ) * (*atoms)[jat].Z;
              jato = jat; 
           }

           double betaij = Bxy[ (*atoms)[iat].atno ][ (*atoms)[jat].atno ];   // resonanace integral for coupling between different atoms
           double Sij    = overlap( bfns[i-1], bfns[j-1] );                   // overlap
           double IPij   = bfns[i-1].ip+bfns[j-1].ip;

           F0(i,j) = betaij*IPij*Sij;

        }
     }
  }

  F0_initialized = TRUE;

}


// Average occupation density matrix
void mindo3::guess_D() {

//  if ( D_initialized ) return;

  D = 0;

  for ( int i=1; i<=nbf; i++ )
      if ((*atoms)[bfns[i-1].natom].atno == 1)
          D(i,i) = (*atoms)[bfns[i-1].natom].Z/1.0;
      else                 
          D(i,i) = (*atoms)[bfns[i-1].natom].Z/4.0;

//  D_initialized = TRUE;

}


// Coulomb-like term for orbitals on the same atom
double mindo3::g(bfn bfi, bfn bfj) {
    int i = bfi.type, j = bfj.type;

    if      ((i==0) && (j==0)) return gss [(*atoms)[bfi.natom].atno];
    else if ((i==0) || (j==0)) return gsp [(*atoms)[bfi.natom].atno];
    else if  (i==j)            return gpp [(*atoms)[bfi.natom].atno];
                               return gppp[(*atoms)[bfi.natom].atno];
}


// Exchange-like term for orbitals on the same atom
double mindo3::h(bfn bfi, bfn bfj) {
    if ((bfi.type==0) || (bfj.type==0)) 
        return hsp [(*atoms)[bfi.natom].atno];
    else       
        return hppp[(*atoms)[bfi.natom].atno];
}              


// One-center corrections to the core fock matrix
void mindo3::calc_F1() {

    F1=0.0;

    for ( int i=1; i<=nbf; i++ ) {
       bfn ibf = bfns[i-1];

       F1(i,i) = 0.5 * g(bfns[i-1],bfns[i-1]) * D(i,i);

       for ( int j=1; j<=nbf; j++ ) { 
          bfn jbf = bfns[j-1];

          if ( (ibf.natom == jbf.natom) && (i!=j) ) {
             double gij = g ( ibf, jbf ),
                    hij = h ( ibf, jbf );
             F1(i,i) += (gij - 0.5*hij) * D(j,j);
             if (i>j) F1(i,j) += 0.5 * (3*hij-gij) * D(i,j);
          }
       }
    }

}


// Two-electron two-center corrections to the core fock matrix
void mindo3::calc_F2() {

    F2=0.0;

    for ( int i=1; i<=nbf; i++ ) {
       int inat = bfns[i-1].natom;

       for ( int j=i+1; j<=nbf; j++ ) {
          int jnat = bfns[j-1].natom;

          if ( inat != jnat ) {
             double gammaij = gamma( (*atoms)[inat], (*atoms)[jnat] );
             F2(i,i) +=         gammaij *   D(j,j);
             F2(j,j) +=         gammaij *   D(i,i);
             F2(i,j) -=  0.25 * gammaij * ( D(i,j) + D(i,j) );
          }
       }
    }

}


// SCF procedure for closed-shell molecules
double mindo3::SCF () { 

    F0.resize(nbf);
    F1.resize(nbf);
    F2.resize(nbf);
    F .resize(nbf);

    calc_F0();

//    guess_D();                                                                // simple guess
    if (!D_initialized) {                                                        // Huckel guess
      ehm::Energy();
      D_initialized = TRUE;
    }

//    SymmetricMatrix D1(nbf); 
//    D1 = D;

    double Eel, Eold = -100 /*,maxdensdiff*/;

//DEBUGSTART("MINDO/3\n");
    for ( SCFit=0; SCFit<MaxSCF; SCFit++ ) {

        calc_F1();
        calc_F2();
        F << F0 + F1 + F2;                                                      // fock matrix F = F0+F1+F2;
//DEBUG("F0:\n");print_mat(F0);
//DEBUG("F0+F2:\n");print_mat(F0+F2);
//DEBUG("F=F0+F1+F2:\n");print_mat(F);
        Eel = 0.5 * trace( D * (F0 + F) );                                      // energy calculation
        Ediff = fabs(Eel-Eold);
        if ( Ediff < TolSCF ) break;                                            // finish SCF if energy difference become small
        Eold = Eel;
        EigenValues( F, Orbe, Orbs );                                           // or Jacobi(F, Dm, U), which is slower, but more reliable
//        D << Orbs.columns(1,nclosed) * Orbs.columns(1,nclosed).t() * 2;         // Form density matrix C*Ct given eigenvectors C[nstart:nstop,:]
        mkdens();
//        maxdensdiff = maximum_absolute_value(D-D1);
//        if ( maxdensdiff < TolSCF ) break;                                      // finish SCF if density matrix deffrence small
//        D1 = D;

    }

//DEBUG("F:\n");print_mat(F);
//DEBUG("D:\n");print_dens();
    return Eel;

}


double mindo3::Energy() {

  double Enuke  = enuke();
  double eref   = refeng();
  double Eel    = SCF();
  double Etot   = Eel+Enuke;
  double Hf     = Etot*ev2kcal+eref;

  return Hf;

}


void mindo3::print_grads(Matrix G) {
  printf("%16s%16s%16s\n","X","Y","Z");
  for( int i=0; i<nat; i++) {
     for ( int j=0;j<3;j++)
        printf("%16.8le",G(i+1,j+1));
     puts("");
  }
}

void mindo3::print_hess(SymmetricMatrix H) {
  for( int i=0; i<nat*3; i++) {
    for( int j=0; j<nat*3; j++) {
        if (i>=j) 
           printf("%13.9lf",H(i+1,j+1)/*/2240.579577*/); 
        else 
           printf("%13s"," ");
     }
     puts("");
  }
}

// Electron density renormalization for SE models:   
// http://www.cachesoftware.com/mopac/Mopac2002manual/node376.html
void mindo3::ElectronDensity_renormalize() {
//  SquareMatrix Orbs_(nbf);
//  Orbs_ << Orbs;
//DEBUGSTART("ElecDen\n");
  SymmetricMatrix S(nbf),Sn(nbf);
  DiagonalMatrix d(nbf);
  SquareMatrix P(nbf);
  int i;

  for ( i=1; i<=nbf; i++) {
    for ( int j=1; j<=nbf; j++) {
      if (i>=j) {
        S(i,j) = overlap ( bfns[i-1], bfns[j-1] );                            // overlap matrix
      }
    }
  }

//DEBUG("S:\n");print_mat(S);
  EigenValues( S, d, P );                                                       // Diagonalize S = P*D*Pt to obtain
  for ( i=1; i<=nbf; i++) d(i,i) = power( d(i,i), -0.5);                    // S^(-1/2) 
  S << P*d*P.t();                                                               // http://seehuhn.de/comp/matrixfn.html
//DEBUG("S=S^(-1/2):\n");print_mat(S);
  Sn << S.i()*S.i();
//DEBUG("S^(-2):\n");print_mat(Sn);
//  S << Sn.i();
//DEBUG("Sn^(-1):\n");print_mat(S);

//DEBUG("Orbs:\n");print_orbs();
  Orbs << S*Orbs*S;                                                             // renormalize MO: Fi = S^(-1/2)*Fi'*S^(-1/2)
//DEBUG("Orbs':\n");print_orbs();

  S.release();
//  Orbs << Orbs_;                                                                // return old orbitals
//print_orbs();

//DEBUG("D:\n");print_dens();
  mkdens();                                                                     // D'
//DEBUG("D':\n");print_dens();

}


// Compute analytic forces on list of atoms
void mindo3::forces(Matrix *forces) {

    int bfi,bfj;
    
    forces->resize(nat,3); *forces=0.0;

    // Loop over all pairs of atoms and compute the force between them
    for (int iat=0, iatb=0; iat<nat; iat++) {
        struct atom_ehm atomi = (*atoms)[iat];
        for (int jat=0, jatb=0; jat<iat; jat++ ) {

            struct atom_ehm atomj = (*atoms)[jat];

            double alpha = axy[atomi.atno][atomj.atno];
            double beta  = Bxy[atomi.atno][atomj.atno];
            double R2 = dist2(atomi.xyz,atomj.xyz);
            double R  = sqrt(R2);
            double c2 = 0.25*power(e2/f03[atomi.atno]+e2/f03[atomj.atno],2);

            for (int dir=0; dir<3; dir++ ){
                double Fij = 0.0;                                               // Force between atoms iat and jat in direction dir

                double delta = atomi.xyz[dir]-atomj.xyz[dir];                   // initialize some constants
                double c1  = delta*atomi.Z*atomj.Z*e2/R;
                double dr1 = e2*delta*power(R2+c2,-1.5);

                if ( ((atomi.atno == 1)                                         // Nuclear repulsion terms
                      && ((atomj.atno == 7) || (atomj.atno == 8)))
                     || ((atomj.atno == 1)
                         && ((atomi.atno == 7) || (atomi.atno == 8))))
                    Fij += -c1*alpha*(1/R2 - R*power(R2+c2,-1.5)                // Special case of NH or OH bonds
                                      + 1/R - 1/sqrt(R2+c2))*exp(-R)
                                      - c1*R*pow(R2+c2,-1.5);
                else
                    Fij += -c1*(1/R2 - R*power(R2+c2,-1.5) + alpha/R 
                                - alpha/sqrt(R2+c2))*exp(-alpha*R)
                                - c1*R*power(R2+c2,-1.5);

                for ( bfi=iatb; bfi<iatb+atomi.nbf; bfi++) {                    // Overlap terms
                    for ( bfj=jatb; bfj<jatb+atomj.nbf; bfj++) {
                        double Dij  = D(bfi+1,bfj+1);
                        double dSij = doverlap( bfns[bfi], bfns[bfj], dir );
                        Fij += 2 * beta * ( bfns[bfi].ip + bfns[bfj].ip ) * Dij * dSij;
                    }
                }
                
                for ( bfj=jatb; bfj<jatb+atomj.nbf; bfj++)                      // Core attraction terms
                    Fij += atomi.Z*D(bfj+1,bfj+1)*dr1;
                for ( bfi=iatb; bfi<iatb+atomi.nbf; bfi++)
                    Fij += atomj.Z*D(bfi+1,bfi+1)*dr1;
                
                for ( bfi=iatb; bfi<iatb+atomi.nbf; bfi++) {                    // Two-electron terms
                   for ( bfj=jatb; bfj<jatb+atomj.nbf; bfj++) {
                        double Dii = D(bfi+1,bfi+1);
                        double Djj = D(bfj+1,bfj+1);
                        double Dij = D(bfi+1,bfj+1);
                        // exchange is the first term, coulomb is second:
                        Fij += 0.5*dr1*power(Dij,2)-dr1*Dii*Djj;
                   }
                }
                
                (*forces)(iat+1,dir+1) += Fij;                                  // Now sum total forces
                (*forces)(jat+1,dir+1) -= Fij;

            }

            jatb += (*atoms)[jat].nbf;
        }

        iatb += (*atoms)[iat].nbf;

    }

    (*forces) *= ev2kcal;                                                       // convert to kcal/mol

}

// Compute numerical forces on list of atoms - Right side finite differencial
// dx = 1.0E-7 ... 1.0E-8 is the best choice according to my experiments
void mindo3::num_forces_right(Matrix *forces, double dx) {

    forces->resize(nat,3); *forces=0.0;

    double fnc1 = Energy();

    for (int iat=0; iat<nat; iat++) {
            for (int dir=0; dir<3; dir++ ){
                
                (*atoms)[iat].xyz[dir] += dx;
                F0_initialized = FALSE;
                double fnc2 = Energy();
                (*atoms)[iat].xyz[dir] -= dx;

                (*forces)(iat+1,dir+1) = (fnc2-fnc1)/dx;
            }
    }
}

// Compute numerical forces on list of atoms - Left side finite differencial
// dx = 1.0E-7 ... 1.0E-8 is the best choice according to my experiments
void mindo3::num_forces_left(Matrix *forces, double dx) {

    forces->resize(nat,3); *forces=0.0;

    double fnc1 = Energy();

    for (int iat=0; iat<nat; iat++) {
            for (int dir=0; dir<3; dir++ ){
                
                (*atoms)[iat].xyz[dir] -= dx;
                F0_initialized = FALSE;
                double fnc2 = Energy();
                (*atoms)[iat].xyz[dir] += dx;

                (*forces)(iat+1,dir+1) = (fnc1-fnc2)/dx;
            }
    }
}

// Compute numerical forces on list of atoms - Central finite differencial
// dx = 1.0E-4 ... 1.0E-8 is the best choice according to my experiments
void mindo3::num_forces_central(Matrix *forces, double dx) {

    forces->resize(nat,3); *forces=0.0;

    for (int iat=0; iat<nat; iat++) {
            for (int dir=0; dir<3; dir++ ){
                
                (*atoms)[iat].xyz[dir] += dx;
                F0_initialized = FALSE;
                double fnc1 = Energy();
                (*atoms)[iat].xyz[dir] -= dx;

                (*atoms)[iat].xyz[dir] -= dx;
                F0_initialized = FALSE;
                double fnc2 = Energy();
                (*atoms)[iat].xyz[dir] += dx;

                (*forces)(iat+1,dir+1) = (fnc1-fnc2)/(dx+dx);
            }
    }
}


// Compute numerical forces on list of atoms - Central four points finite differencial
void mindo3::num_forces_central4(Matrix *forces, double dx) {

    forces->resize(nat,3); *forces=0.0;

    for (int iat=0; iat<nat; iat++) {
            for (int dir=0; dir<3; dir++ ){
                
                (*atoms)[iat].xyz[dir] -= 2*dx;
                F0_initialized = FALSE;
                double fnc1 = Energy();
                (*atoms)[iat].xyz[dir] += 2*dx;

                (*atoms)[iat].xyz[dir] -= dx;
                F0_initialized = FALSE;
                double fnc2 = Energy();
                (*atoms)[iat].xyz[dir] += dx;

                (*atoms)[iat].xyz[dir] += dx;
                F0_initialized = FALSE;
                double fnc3 = Energy();
                (*atoms)[iat].xyz[dir] -= dx;

                (*atoms)[iat].xyz[dir] += 2*dx;
                F0_initialized = FALSE;
                double fnc4 = Energy();
                (*atoms)[iat].xyz[dir] -= 2*dx;

                (*forces)(iat+1,dir+1) = (fnc1-8*fnc2+8*fnc3-fnc4)/(12*dx);
//                (*forces)(iat+1,dir+1) = (-3*fnc1+8*fnc2-24*fnc3+19*fnc4)/(12*dx);

            }
    }
}


// Compute numerical hessian on list of atoms - Right side finite differencial
void mindo3::num_hessian_right(SymmetricMatrix *hess, double dx) {

    Matrix fr1, fr2;

    int i, j, iat;
    hess->resize(nat*3); *hess=0.0;

    Energy();
    forces( &fr1 );


    DiagonalMatrix xi(nat*3);
    for (iat=0; iat<nat; iat++) {
        xi(iat*3+1) = -fr1(iat+1,1);
        xi(iat*3+2) = -fr1(iat+1,2);
        xi(iat*3+3) = -fr1(iat+1,3);
    }                       

    for (i=0; i<nat*3; i++) {
        (*hess)(i+1,i+1) = 1.0;
    }
    
    for (int k=0; k<1; k++) {
    
        for (iat=0; iat<1/*nat*/; iat++) {
            (*atoms)[iat].xyz[0] += dx;
            (*atoms)[iat].xyz[1] += dx;  
            (*atoms)[iat].xyz[2] += dx;  
        }

        Energy();
        forces(&fr2);

        DiagonalMatrix dg(nat*3);
        for (i=0; i<nat*3; i++) {
            int iat=i/3, dir=i%3;
            dg(i+1) = fr2(iat+1,dir+1) - fr1(iat+1,dir+1);
        }

        DiagonalMatrix hdg(nat*3); hdg = 0.0;
        for (i=1; i<=nat*3; i++) {
           for (j=1; j<=nat*3; j++) {
              hdg(i) += (*hess)(i,j)*dg(j);
           }
        }

        double fac = 0.0;
        double fae = 0.0;

        for (i=1; i<=nat*3; i++) {
           fac += dg(i)*xi(i);
           fae += dg(i)*hdg(i);
        }

        fac = 1/fac;
        double fad = 1/fae;

        for (i=1; i<=nat*3; i++) {
           dg(i) = fac*xi(i)-fad*hdg(i);
        }

        for (i=1; i<=nat*3; i++) {
          for (j=1; j<=nat*3; j++) {
             (*hess)(i,j) += fac*xi(i)*xi(j) - fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j);
          }
        }
        

        for (int iat=0; iat<1/*nat*/; iat++) {
            (*atoms)[iat].xyz[0] -= dx;
            (*atoms)[iat].xyz[1] -= dx;
            (*atoms)[iat].xyz[2] -= dx;
        }

    }

}


// Compute numerical hessian on list of atoms - Both side finite differencial
void mindo3::num_hessian_central(SymmetricMatrix *hess, double dx) {

    Matrix fr1, fr2;

    int i, j;
    hess->resize(nat*3); *hess=0.0;

    for (i=0; i<nat*3; i++) {
        int iat=i/3, dir=i%3;
        
        (*atoms)[iat].xyz[dir] += dx;
        Energy();
        forces(&fr2);
        (*atoms)[iat].xyz[dir] -= dx;

        (*atoms)[iat].xyz[dir] -= dx;
        Energy();
        forces(&fr1);
        (*atoms)[iat].xyz[dir] += dx;

        for (j=0; j<=i; j++) {
           int jat=j/3, djr=j%3;
           (*hess)(i+1,j+1) = (fr2(jat+1,djr+1)-fr1(jat+1,djr+1))/(dx+dx);
        }
    }
}
