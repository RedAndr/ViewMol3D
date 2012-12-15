
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#ifdef WIN32
#include <windows.h>
#endif

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstdio>

#include <vector>
#include <string>


using namespace std;

#include "newmatap.h"                   // newmat headers including advanced functions
#include "atom_rads.h"
#include "readfile.h"

#define len(x,y,z,xx,yy,zz) sqrt(sqr((x)-(xx))+sqr((y)-(yy))+sqr((z)-(zz)))
#define sqr(a)  ((a)*(a))

ReadMolecule::ReadMolecule() :
    tstrn    (0),
    ircsave  (0),            
    freqss   (0),             
    gradss   (0),
    NumData  (0),
    CurrentPoint (0),
    ObjCount (0),             
    cread    (0),
    bread    (0)
//    atoms(NULL), bonds(NULL), grads(NULL), freqs(NULL)
{
    yes[0]  = 0;
    yes[1]  = 0;
    yes[2]  = 0;
    yes[3]  = 0;
}

int ReadMolecule::OutFile (char* fname, int strn) {
  char s[128],stype[128],ss[128];
  int i,a[6],ii,j,in[3],stdor;
  int strcur=0;
  float it,x,y,z,b[6];
  FILE *f, *firc;
  struct atom atom;
  struct freq freq;
  struct grad grad;

  bread=0;

  f=fopen(fname,"r");
  if(f==NULL) return -1;
//if(debug)printf("Open GAMESS/GAUSSIAN98 OUT-file: %s\n",fname);

  if(ircsave) firc=fopen("irc.out","w");

  do {
    fgets(s,100,f);
    if ((strstr(s,"ATOM      ATOMIC                      COORDINATES (BOHR)")!=NULL)) {      // GAMESS coordinates
      fgets(s,120,f);
      for (int i=0;;i++) {
        fgets(s,120,f);
        if (s[0]==0x0A) break;
        sscanf(s,"%s %f %f %f %f\n",&stype,&it,&x,&y,&z);
        atom.x=(double)x*0.529177249;  
        atom.y=(double)y*0.529177249;  
        atom.z=(double)z*0.529177249;
        atom.type=(int)it;
        if ( cread == 1 ) atoms[i] = atom;
        else atoms.push_back( atom );
      }
      cread=1;
      strcur++; 
      if ( strcur==strn ) { 
        fclose(f); 
        return 0; 
      }
    };

    if (strstr(s,"COORDINATES OF ALL ATOMS ARE (ANGS)")!=NULL) {      // GAMESS coordinates
      fgets(s,120,f); fgets(s,120,f);
      for (int i=0;;i++) {
        fgets(s,120,f);
        if ((s[0]==0x0A)||(strstr(s," COORDINATES OF")!=NULL)) break;
        sscanf(s,"%s %f %f %f %f\n",&stype,&it,&x,&y,&z);
        atom.type=(int)it;
        atom.x=(double)x;  
        atom.y=(double)y;  
        atom.z=(double)z;
        if ( cread == 1 ) atoms[i] = atom;
        else atoms.push_back( atom );
      }
      cread = 1;
      strcur++; 
      if(strcur==strn) { 
        fclose(f); 
        return 0; 
      }
    };

/*    if (strstr(s,"NET REACTION COORDINATE UP TO THIS POINT")!=NULL) {
      sscanf(s,"%s %s %s %s %s %s %s %s %f\n",&stype,&stype,&stype,&stype,&stype,&stype,&stype,&stype,&x);
      CurrentIRC=x;
if(debug)printf("Current IRC=%f\n",CurrentIRC);
      if(ircsave) {
        fprintf(firc,"NET REACTION COORDINATE UP TO THIS POINT = %7.5f\n",CurrentIRC);      
      }
      int flag=1;
      for(i=0;i<NumData;i++) if(CurrentIRC==GraphData[i][0]) flag=0;
      if(flag){
        GraphData[NumData][0]=CurrentIRC;
        GraphData[NumData][1]=CurrentEnergy;
        NumData++;
      };
    }*/

    if (strstr(s,"SCF Done:")!=NULL) {
          sscanf(strstr(s,"="),"%s %lf",&stype,&CurrentEnergy);
//if(debug)printf("Current energy=%15.9lf\n",CurrentEnergy);
      if(ircsave) fprintf(firc,"SCF Done:  E = %15.9lf A.U.\n",CurrentEnergy);      
    }

    if (strstr(s,"Item               Value     Threshold")!=NULL) {
      fgets(s,120,f); if (strstr(s,"YES")!=NULL) yes[0]=1; else yes[0]=0;
      fgets(s,120,f); if (strstr(s,"YES")!=NULL) yes[1]=1; else yes[1]=0;
      fgets(s,120,f); if (strstr(s,"YES")!=NULL) yes[2]=1; else yes[2]=0;
      fgets(s,120,f); if (strstr(s,"YES")!=NULL) yes[3]=1; else yes[3]=0;
      if(yes[0]&&yes[1]&&yes[2]&&yes[3]&&ircsave) fprintf(firc,"Item               Value     Threshold\nYES\nYES\nYES\nYES\n");      
//if(debug)printf("%d %d %d %d\n",yes[0],yes[1],yes[2],yes[3]);

/*      if(yes[0]&&yes[1]&&yes[2]&&yes[3]){
if(debug)printf("Found minimum #%3u Energy:%f\n",CurrentPoint,CurrentEnergy);
          GraphData[NumData][0]=(double)strcur;//CurrentPoint;
          GraphData[NumData][1]=CurrentEnergy;
          CurrentPoint++;
          NumData++;
      }*/
    }

    if (strstr(s,"R6    R(3,8)")!=NULL && CurrentEnergy!=0) {
      sscanf(s,"%s %s %s %lf\n",&stype,&stype,&stype,&CurrentIRC);
//if(debug)printf("Current IRC=%5.3f CurrentEnergy=%15.9f\n",CurrentIRC,CurrentEnergy);
      if(ircsave) fprintf(firc,"R6    R(3,8) %5.3lf\n",CurrentIRC);      
      int flag=1;
      for(i=0;i<NumData;i++) if(CurrentIRC==GraphData[i][0]) flag=0;
      if(flag){
        GraphData[NumData][0]=CurrentIRC;
        GraphData[NumData][1]=CurrentEnergy;
        NumData++;
      }
    }

    if (/*strstr(s,"Standard orientation:")!=NULL *||*/
        strstr(s,"Input orientation:"   )!=NULL /*||
        strstr(s,"Z-Matrix orientation:")!=NULL*/
       ) {

      if(strstr(s,"Standard orientation:")!=NULL) stdor=1; else stdor=0;

//if(debug)puts("Found G98/03 cartesian coordinates, reading ...");
      fgets(s,120,f); fgets(s,120,f); fgets(s,120,f); fgets(s,120,f);
      
      for (int i=0;;i++) {
        fgets(s,120,f);
        if ((s[0]==0x0A)||(strstr(s,"------------------")!=NULL)) break;
        sscanf(s,"%d %d %d %f %f %f\n",&in[0],&in[1],&in[2],&x,&y,&z);
        atom.x=(double)x;  
        atom.y=(double)y;
        atom.z=(double)z;
        atom.type=in[1];
        if ( cread == 1 ) atoms[i] = atom;
        else atoms.push_back( atom );
      }

      cread=1;
      strcur++; 

      if(strcur==strn) { 
//        if(bread==0){ CreateBond(); bread=1; };
        fclose(f);
        return 0; 
      }

      if(ircsave && stdor && yes[0] && yes[1] && yes[2] && yes[3]) {
         fputs("                         Standard orientation:\n",firc);
         fputs(" ---------------------------------------------------------------------\n",firc);
         fputs(" Center     Atomic     Atomic              Coordinates (Angstroms)\n",firc);
         fputs(" Number     Number      Type              X           Y           Z\n",firc);
         fputs(" ---------------------------------------------------------------------\n",firc);
         for(i=0;i<atoms.size();i++)
           fprintf(firc,"%3d %3d %3d %9.6f %9.6f %9.6f\n",i+1,atoms[i].type,i,atoms[i].x,atoms[i].y,atoms[i].z);
         fputs(" ---------------------------------------------------------------------\n",firc);
         fputs("\n\n",firc);
      }

    };

/*    if (strstr(s,"BOND ORDER")!=NULL) {
if(debug)puts("Found bond orders, reading ...");
      fgets(s,120,f); fgets(s,120,f); fgets(s,120,f); fgets(s,120,f);
      i=0;
      for (;;) {
        fgets(s,120,f);
if(debug)puts(s);
        if (s[0]==0x0A) break;
        a[0]=0; a[2]=0; a[4]=0;
        sscanf(s,"%d %d %f %f %d %d %f %f %d %d %f %f \n",
                                 &a[0],&a[1],&b[0],&b[1],
                                 &a[2],&a[3],&b[2],&b[3],
                                 &a[4],&a[5],&b[4],&b[5]);
        if (bread==0) bonds = (struct bond *) realloc (bonds, (i+3)*sizeof(struct bond));
        if((a[0]!=0)&&(b[1]>=BondTresh)) {
          bonds[i].a=a[0]; bonds[i].b=a[1]; bonds[i].o=b[1]/10;
if(debug)printf("bond:%2d %3d - %3d : %5.3f\n",i,bonds[i].a,bonds[i].b,b[1]); i++;}
        if((a[2]!=0)&&(b[3]>=BondTresh)) {
          bonds[i].a=a[2]; bonds[i].b=a[3]; bonds[i].o=b[3]/10;
if(debug)printf("bond:%2d %3d - %3d : %5.3f\n",i,bonds[i].a,bonds[i].b,b[3]); i++;}
        if((a[4]!=0)&&(b[5]>=BondTresh)) {
          bonds[i].a=a[4]; bonds[i].b=a[5]; bonds[i].o=b[5]/10;
if(debug)printf("bond:%2d %3d - %3d : %5.3f\n",i,bonds[i].a,bonds[i].b,b[5]); i++;}
      }
      NumBonds=i;
      if(debug)printf("NumBonds=%d\n",NumBonds);
      bread=1;
    }*/

    if (strstr(s,"TOTAL MULLIKEN")!=NULL) {
//if(debug)puts("Found mulliken charges, reading ...");
      fgets(s,120,f);
      i=0;
      for (;;) {
        fgets(s,120,f);
        if (s[0]==0x0A) break;
        sscanf(s,"%d %s %f %f %f %f\n",&a[0],stype,&b[0],&b[1],&b[2],&b[3]);
        atoms[i].c=b[1]/2;
//        if(debug)printf("%3d %7.4f\n",i,b[1]);
        i++;
      }
    }

    if (strstr(s,"MOPAC CHARGES")!=NULL) {
//if(debug)puts("Found MOPAC charges, reading ...");
      fgets(s,120,f);  fgets(s,120,f);
      i=0;
      for (;;) {
        fgets(s,120,f);
        if (s[0]==0x0A) break;
        sscanf(s,"%d %s %f %f\n",&a[0],stype,&b[0],&b[1]);
        atoms[i].c=b[0]/2;
//        if(debug)printf("%3d %7.4f\n",i,b[1]);
        i++;
      }
    }

    if (strstr(s,"Total atomic charges:")!=NULL) {
//if(debug)puts("Found G98 charges, reading ...");
      fgets(s,120,f);
      for (i=0;;i++) {
        fgets(s,120,f);
        if ((s[0]==0x0A)||(strstr(s,"Sum of Mulliken charges")!=NULL)) break;
        sscanf(s,"%d %s %f\n",&a[0],&stype,&b[0]);
        atoms[i].c=b[0]/2;
//if(debug)printf("%3d %7.4f\n",i,b[0]);
      }
    }

    if (strstr(s,"FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2")!=NULL) {
//if(debug)puts("Found GAMESS frequencies, reading ...");

      struct freq freqa[5];
      struct atom atoma[5];

      for(i=1;i<=atoms.size()*3-6;i+=5) {
        fgets(s,90,f);          // Space
        fgets(s,90,f);          // Numbers of Freqs
        if (strstr(s,"***************")!=NULL) for(j=1;j<=5;j++) fgets(s,90,f);
        fgets(s,90,f);          // Freqs
        ss[8]=0;
        sscanf(strncpy(ss,s+22,8),"%f",&it); freqa[0].freqv=(double)it;
        sscanf(strncpy(ss,s+34,8),"%f",&it); freqa[1].freqv=(double)it;
        sscanf(strncpy(ss,s+46,8),"%f",&it); freqa[2].freqv=(double)it;
        sscanf(strncpy(ss,s+58,8),"%f",&it); freqa[3].freqv=(double)it;
        sscanf(strncpy(ss,s+70,8),"%f",&it); freqa[4].freqv=(double)it;
        ss[1]=0;
        if (strstr(strncpy(ss,s+31,1),"I")!=NULL) freqa[0].img=1; else  freqa[0].img=0;
        if (strstr(strncpy(ss,s+43,1),"I")!=NULL) freqa[1].img=1; else  freqa[1].img=0;
        if (strstr(strncpy(ss,s+55,1),"I")!=NULL) freqa[2].img=1; else  freqa[2].img=0;
        if (strstr(strncpy(ss,s+67,1),"I")!=NULL) freqa[3].img=1; else  freqa[3].img=0;
        if (strstr(strncpy(ss,s+79,1),"I")!=NULL) freqa[4].img=1; else  freqa[4].img=0;
        fgets(s,90,f);          // Intens
        if (strstr(s,"INTENSITY:")!=NULL) {
          sscanf(strncpy(ss,s+22,8),"%f",&it); freqa[0].inten=(double)it;
          sscanf(strncpy(ss,s+34,8),"%f",&it); freqa[1].inten=(double)it;
          sscanf(strncpy(ss,s+46,8),"%f",&it); freqa[2].inten=(double)it;
          sscanf(strncpy(ss,s+58,8),"%f",&it); freqa[3].inten=(double)it;
          sscanf(strncpy(ss,s+70,8),"%f",&it); freqa[4].inten=(double)it;
          fgets(s,90,f);          // Space
        } else for(j=0;j<5;j++) freqa[j].inten=0;
        ss[11]=0;

        freqa[0].coords.erase( freqa[0].coords.begin(), freqa[0].coords.end() );
        freqa[1].coords.erase( freqa[1].coords.begin(), freqa[1].coords.end() );
        freqa[2].coords.erase( freqa[2].coords.begin(), freqa[2].coords.end() );
        freqa[3].coords.erase( freqa[3].coords.begin(), freqa[3].coords.end() );
        freqa[4].coords.erase( freqa[4].coords.begin(), freqa[4].coords.end() );
        
        for(j=1;j<=atoms.size();j++) {
          fgets(s,90,f);        // x
          sscanf(strncpy(ss,s+21,11),"%f",&it); atoma[0].x=(double)it;
          sscanf(strncpy(ss,s+33,11),"%f",&it); atoma[1].x=(double)it;
          sscanf(strncpy(ss,s+45,11),"%f",&it); atoma[2].x=(double)it;
          sscanf(strncpy(ss,s+57,11),"%f",&it); atoma[3].x=(double)it;
          sscanf(strncpy(ss,s+69,11),"%f",&it); atoma[4].x=(double)it;
          fgets(s,90,f);        // y
          sscanf(strncpy(ss,s+21,11),"%f",&it); atoma[0].y=(double)it;
          sscanf(strncpy(ss,s+33,11),"%f",&it); atoma[1].y=(double)it;
          sscanf(strncpy(ss,s+45,11),"%f",&it); atoma[2].y=(double)it;
          sscanf(strncpy(ss,s+57,11),"%f",&it); atoma[3].y=(double)it;
          sscanf(strncpy(ss,s+69,11),"%f",&it); atoma[4].y=(double)it;
          fgets(s,90,f);        // z
          sscanf(strncpy(ss,s+21,11),"%f",&it); atoma[0].z=(double)it;
          sscanf(strncpy(ss,s+33,11),"%f",&it); atoma[1].z=(double)it;
          sscanf(strncpy(ss,s+45,11),"%f",&it); atoma[2].z=(double)it;
          sscanf(strncpy(ss,s+57,11),"%f",&it); atoma[3].z=(double)it;
          sscanf(strncpy(ss,s+69,11),"%f",&it); atoma[4].z=(double)it;
          freqa[0].coords.push_back( atoma[0] );
          freqa[1].coords.push_back( atoma[1] );
          freqa[2].coords.push_back( atoma[2] );
          freqa[3].coords.push_back( atoma[3] );
          freqa[4].coords.push_back( atoma[4] );
        }
        for(j=1;j<=10;j++) fgets(s,90,f);          // Garbage

        freqs.push_back( freqa[0] );
        freqs.push_back( freqa[1] );
        freqs.push_back( freqa[2] );
        freqs.push_back( freqa[3] );
        freqs.push_back( freqa[4] );
      }
      freqss=1;
    }
    
    if (strstr(s,"reduced masses (AMU), force constants (mDyne/A) and normal coordinates:")!=NULL) {
//if(debug)puts("Found G98 frequencies, reading ...");

      struct freq freq1;
      struct freq freq2;
      struct freq freq3;
      struct atom atom;

      for(i=1;i<=atoms.size()*3-6;i+=3) {
        fgets(s,90,f);          // Numbers of Freqs
        fgets(s,90,f);          // A'
        fgets(s,90,f);          // Freqs
        ss[9]=0;
        sscanf(strncpy(ss,s+16,10),"%f",&it); if(it<0) {it=-it; freq1.img=1;} else freq1.img=0; freq1.freqv=(double)it;
        sscanf(strncpy(ss,s+39,10),"%f",&it); if(it<0) {it=-it; freq2.img=1;} else freq2.img=0; freq2.freqv=(double)it;
        sscanf(strncpy(ss,s+62,10),"%f",&it); if(it<0) {it=-it; freq3.img=1;} else freq3.img=0; freq3.freqv=(double)it;
        fgets(s,90,f);          // Red. masses
        fgets(s,90,f);          // Frc consts
        fgets(s,90,f);          // IR Inten
        sscanf(strncpy(ss,s+16,10),"%f",&it); freq1.inten=(double)it;
        sscanf(strncpy(ss,s+39,10),"%f",&it); freq2.inten=(double)it;
        sscanf(strncpy(ss,s+62,10),"%f",&it); freq3.inten=(double)it;
        fgets(s,90,f);          // Raman Activ
        fgets(s,90,f);          // Depolar
        fgets(s,90,f);          // Atom AN
        ss[5]=0;
        freq1.coords.erase( freq1.coords.begin(), freq1.coords.end() );
        freq2.coords.erase( freq2.coords.begin(), freq2.coords.end() );
        freq3.coords.erase( freq3.coords.begin(), freq3.coords.end() );
        for(j=1;j<=atoms.size();j++) {
          fgets(s,90,f);        // x
          sscanf(strncpy(ss,s+12,5),"%f",&it); atom.x=(double)it;
          sscanf(strncpy(ss,s+19,5),"%f",&it); atom.y=(double)it;
          sscanf(strncpy(ss,s+26,5),"%f",&it); atom.z=(double)it;
          freq1.coords.push_back( atom );
          sscanf(strncpy(ss,s+35,5),"%f",&it); atom.x=(double)it;
          sscanf(strncpy(ss,s+42,5),"%f",&it); atom.y=(double)it;
          sscanf(strncpy(ss,s+49,5),"%f",&it); atom.z=(double)it;
          freq2.coords.push_back( atom );
          sscanf(strncpy(ss,s+58,5),"%f",&it); atom.x=(double)it;
          sscanf(strncpy(ss,s+65,5),"%f",&it); atom.y=(double)it;
          sscanf(strncpy(ss,s+72,5),"%f",&it); atom.z=(double)it;
          freq3.coords.push_back( atom );
        }
        freqs.push_back( freq1 );
        freqs.push_back( freq2 );
        freqs.push_back( freq3 );
      }
      freqss=1;

    }

    if (strstr(s,"and normal coordinates:")!=NULL) {
//if(debug)puts("Found G03 frequencies, reading ...");

      struct freq freq1;
      struct freq freq2;
      struct freq freq3;
      struct atom atom;

      for(i=1;i<=atoms.size()*3-6;i+=3) {
        fgets(s,90,f);          // Numbers of Freqs
        fgets(s,90,f);          // A'
        fgets(s,90,f);          // Freqs
        ss[9]=0;
        sscanf(strncpy(ss,s+16,10),"%f",&it); if(it<0) {it=-it; freq1.img=1;} else freq1.img=0; freq1.freqv=(double)it;
        sscanf(strncpy(ss,s+39,10),"%f",&it); if(it<0) {it=-it; freq2.img=1;} else freq2.img=0; freq2.freqv=(double)it;
        sscanf(strncpy(ss,s+62,10),"%f",&it); if(it<0) {it=-it; freq3.img=1;} else freq3.img=0; freq3.freqv=(double)it;
        fgets(s,90,f);          // Red. masses                                                            
        fgets(s,90,f);          // Frc consts                                                             
        fgets(s,90,f);          // IR Inten
        sscanf(strncpy(ss,s+16,10),"%f",&it); freq1.inten=(double)it;
        sscanf(strncpy(ss,s+39,10),"%f",&it); freq2.inten=(double)it;
        sscanf(strncpy(ss,s+62,10),"%f",&it); freq3.inten=(double)it;
        fgets(s,90,f);          // Atom AN
        ss[5]=0;
        freq1.coords.erase( freq1.coords.begin(), freq1.coords.end() );
        freq2.coords.erase( freq2.coords.begin(), freq2.coords.end() );
        freq3.coords.erase( freq3.coords.begin(), freq3.coords.end() );
        for(j=1;j<=atoms.size();j++) {
          fgets(s,90,f);        // x
          sscanf(strncpy(ss,s+12,5),"%f",&it); atom.x=(double)it;
          sscanf(strncpy(ss,s+19,5),"%f",&it); atom.y=(double)it;
          sscanf(strncpy(ss,s+26,5),"%f",&it); atom.z=(double)it;
          freq1.coords.push_back( atom );
          sscanf(strncpy(ss,s+35,5),"%f",&it); atom.x=(double)it;
          sscanf(strncpy(ss,s+42,5),"%f",&it); atom.y=(double)it;
          sscanf(strncpy(ss,s+49,5),"%f",&it); atom.z=(double)it;
          freq2.coords.push_back( atom );
          sscanf(strncpy(ss,s+58,5),"%f",&it); atom.x=(double)it;
          sscanf(strncpy(ss,s+65,5),"%f",&it); atom.y=(double)it;
          sscanf(strncpy(ss,s+72,5),"%f",&it); atom.z=(double)it;
          freq3.coords.push_back( atom );
        }
        freqs.push_back( freq1 );
        freqs.push_back( freq2 );
        freqs.push_back( freq3 );
      }
      freqss=1;
    }
    
/*    if (strstr(s,"Eigenvectors required to have negative eigenvalues:")!=NULL){

      fgets(s,120,f);
      if (strstr(s,"X1        Y1        Z1")!=NULL){
        struct freq freq1;
   
        freq1.img=1;
        for(i=1;i<=atoms.size()*3/5+1;i++) {
          if(i!=1) fgets(s ,90,f);          // XnYnZn
          fgets(s2,90,f);                   // Freqs
          ss[2]=0; sscanf(strncpy(ss,s +27,2),"%d",&ii); // Atom number
          ss[8]=0; sscanf(strncpy(ss,s2+23,8),"%f",&it); // Vector
          if(s[26]==88) {freq1.coords[ii].x=(double)it;}
          if(s[26]==89) {freq1.coords[ii].y=(double)it;}
          if(s[26]==90) {freq1.coords[ii].z=(double)it;}

          ss[2]=0; sscanf(strncpy(ss,s +37,2),"%d",&ii); // Atom number
          ss[8]=0; sscanf(strncpy(ss,s2+33,8),"%f",&it); // Vector
          if(s[36]==88) {freq1.coords[ii].x=(double)it;}
          if(s[36]==89) {freq1.coords[ii].y=(double)it;}
          if(s[36]==90) {freq1.coords[ii].z=(double)it;}

          ss[2]=0; sscanf(strncpy(ss,s +47,2),"%d",&ii); // Atom number
          ss[8]=0; sscanf(strncpy(ss,s2+43,8),"%f",&it); // Vector
          if(s[46]==88) {freq1.coords[ii].x=(double)it;}
          if(s[46]==89) {freq1.coords[ii].y=(double)it;}
          if(s[46]==90) {freq1.coords[ii].z=(double)it;}

          ss[2]=0; sscanf(strncpy(ss,s +57,2),"%d",&ii); // Atom number
          ss[8]=0; sscanf(strncpy(ss,s2+53,8),"%f",&it); // Vector
          if(s[56]==88) {freq1.coords[ii].x=(double)it;}
          if(s[56]==89) {freq1.coords[ii].y=(double)it;}
          if(s[56]==90) {freq1.coords[ii].z=(double)it;}

          ss[2]=0; sscanf(strncpy(ss,s +67,2),"%d",&ii); // Atom number
          ss[8]=0; sscanf(strncpy(ss,s2+63,8),"%f",&it); // Vector
          if(s[66]==88) {freq1.coords[ii].x=(double)it;}
          if(s[66]==89) {freq1.coords[ii].y=(double)it;}
          if(s[66]==90) {freq1.coords[ii].z=(double)it;}
        }

//        freqs.push_back( freq1 );
          freqs[0]=freq1;

      }
      freqss=1;
    }*/
    
    if (strstr(s,"GRADIENT (HARTREE/BOHR)")!=NULL) {
//if(debug)puts("Found GAMESS gradients, reading ...");
      fgets(s,90,f);          // Space
      fgets(s,90,f);          // Space
      fgets(s,90,f);          // Space
      for(i=0;i<atoms.size();i++) {
        fgets(s,90,f);          // Grad
        sscanf(s,"%d %s %f %f %f %f\n",&ii,&stype,&it,&x,&y,&z);
        grad.znuc=int(it); 
        grad.x=x; 
        grad.y=y; 
        grad.z=z; 
        if (gradss == 1) grads[i] = grad;
        else grads.push_back( grad );
      }
      gradss=1;
    }
    
    if (strstr(s,"Forces (Hartrees/Bohr)")!=NULL) {
//if(debug)puts("Found Gaussian gradients, reading ...");
      fgets(s,90,f);          // Space
      fgets(s,90,f);          // Space
      for(i=0;i<atoms.size();i++) {
        fgets(s,90,f);          // Grad
        sscanf(s,"%d %d %f %f %f\n",&ii,&j,&x,&y,&z);
        grad.znuc=j; 
        grad.x=x; 
        grad.y=y; 
        grad.z=z; 
        if (gradss == 1) grads[i] = grad;
        else grads.push_back( grad );
      }
      gradss=1;
    }
    
    if (strstr(s,"CARTESIAN COORDINATES")!=NULL) {
//if(debug)puts("Found MOPAC cartesian coordinates, reading ...");
      fgets(s,120,f); fgets(s,120,f); fgets(s,120,f);
      i=0;
      for (;;) {
        fgets(s,120,f);
        if (s[0]==0x0A) break;
        sscanf(s,"%u %s %f %f %f\n",&ii,&s,&x,&y,&z);
        atom.x=(double)x;  
        atom.y=(double)y;  
        atom.z=(double)z;
                       atom.type= 1;
        if (s[0]=='H') atom.type= 1;
        if (s[0]=='B') atom.type= 5;
        if (s[0]=='C') atom.type= 6;
        if (s[0]=='N') atom.type= 7;
        if (s[0]=='O') atom.type= 8;
        if (s[0]=='F') atom.type= 9;
        if (s[0]=='P') atom.type=15;
        if (s[0]=='S') atom.type=16;
        if (s[0]=='I') atom.type=53;
        if((s[0]=='S') && (s[1]=='i')) 
                       atom.type=14;
        if((s[0]=='C') && (s[1]=='l')) 
                       atom.type=17;
        atoms.push_back( atom );
        if(atom.type==-1) break;
        i++;
      }
      cread=1;
      strcur++; 
      if(strcur==strn) { 
        fclose(f); 
        return 0; 
      }
    }
/*
    if (strstr(s,"Condensed to atoms (all electrons):")!=NULL) {
      if(pop==NULL) pop = (double *)calloc(Numatoms*Numatoms,sizeof(double));
if(debug)puts("Found G98 SCF population matrix...");
      for(i=0;i<(Numatoms/6+1);i++) {
        fgets(s,120,f);   // Numbers
        for(j=0;j<Numatoms;j++) {
          fgets(s,120,f);
          sscanf(s,"%d %s  %f %f %f %f %f %f\n",&ii,&s,&b[0],&b[1],&b[2],&b[3],&b[4],&b[5]);
          for(k=0;k<(min(Numatoms-i*6,6));k++) {
            pop[j*Numatoms+i*6+k]=b[k];
          }
if(debug)printf("%d %d : %d `%s'  %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n",i,j,ii,s,b[0],b[1],b[2],b[3],b[4],b[5]);
        }
      }
if(debug)for(i=0;i<Numatoms;i++) { for(j=0;j<Numatoms;j++) printf("%10.6f",pop[i*Numatoms+j]); printf("\n"); }
    }
*/
    if (strstr(s,"VM3_Graph")!=NULL) {
      sscanf(s,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&stype,
                                                    &ObjData[ObjCount][0],                           // Type: 1: sphere, 2: cylinder
        &ObjData[ObjCount][1],&ObjData[ObjCount][2],&ObjData[ObjCount][3],                           // x y z  coordinates
        &ObjData[ObjCount][4],&ObjData[ObjCount][5],&ObjData[ObjCount][6],&ObjData[ObjCount][7],     // x y z  angles
        &ObjData[ObjCount][8],&ObjData[ObjCount][9],&ObjData[ObjCount][10],                          // parameters: 1:rad,sl,st 
        &ObjData[ObjCount][11],&ObjData[ObjCount][12],                                               // 2:rad1,rad2,height,sl,st
        &ObjData[ObjCount][13],&ObjData[ObjCount][14],&ObjData[ObjCount][15],&ObjData[ObjCount][16]);// r g b alpha  colors
      ObjCount++;
//printf("Read VM3_Graph %f\n",ObjData[ObjCount][0]);
//if(fabs(ObjData[ObjCount][0]-1)<0.001) 
    }
/*
    if (strstr(s,"NR  ATOM    CHARGE       X              Y              Z")!=NULL) {
if(debug)printf("Found MolPro XYZ coordinates in Bohr");
      fgets(s,80,f);
      for (;;) {
        fgets(s,120,f);
        if (s[0]==0x0A) break;
        sscanf(s,"%u %s %f %f %f %f\n",&ii,&s,&it,&x,&y,&z);
        atom.x=(double)x;  
        atom.y=(double)y;  
        atom.z=(double)z;
        atom.type=(double)it;
        atoms.push_back( atom );
        i++;
      }
      cread=1;
    }
*/    
    if (strstr(s,"Dump information in style XYZ to output")!=NULL) {            // Found MolPro XYZ coordinates
      fgets(s,80,f); 
      fgets(s,80,f); 
      fgets(s,80,f); 

      char an[8];
      int numatoms;

      sscanf(s,"%d\n",&numatoms);
      fgets(s,80,f); //if(debug)puts(s);
//char zZz[256]; DEBUGSTART("ReadMolPro\n");
      for(i=0;i<numatoms;i++) {
        s[0]=' '; fgets(s+1,80,f);
        sscanf(s,"%s %f %f %f", &s, &x, &y, &z);
        atom.x=(double)x;
        atom.y=(double)y;
        atom.z=(double)z;
        atom.type = 1;                                                           // default atom

        if(atoi(s)>0) 
           atoms[i].type=atoi(s); 
        else
           for( j=0; j<at_namc; j++ ) {

              strcpy(an,at_nam[j]);
              strupr(an);

              if (s[0]==an[0])
                  if (s[1]==an[1]) { 
                    atom.type=j+1; break;
                  } else
                      if (an[1]==' ')
                         atom.type=j+1;
           }
        atoms.push_back( atom );
      }
      cread=1;
    }

  } while (feof(f)==0);

  fclose(f);
  if (ircsave) fclose(firc);

  if (cread==0){ /*puts("Cartesian coordinates not found.");*/ return -1;}

  tstrn=strcur;

  return 0;
}


int ReadMolecule::ALCHEMY (char* fname) {

  char s[80];
  int i,ii,j;
  FILE *f;
  double fl;
  float x,y,z;
  struct atom atom;
  struct bond bond;

  f=fopen(fname,"r");  if(f==NULL) return -1;
  fgets(s,80,f);
  int numatoms, numbonds;
  sscanf(s,"%d %s %d %s %d %s %s\n",&numatoms,&s,&numbonds,&s,&ii,&s,&s);

  for( i=0; i<numatoms; i++) {
    fgets(s,80,f);
    sscanf(s,"%d %s %f %f %f %f\n",
              &ii, &s, &x, &y, &z, &fl);
    atom.x=(double)x;
    atom.y=(double)y;
    atom.z=(double)z;

    atom.type=-1;
    for(j=0;j<at_namc;j++)
      if (s[0]==at_nam[j][0])
        if (s[1]==at_nam[j][1]) {atom.type=j+1; break; } else
          if (at_nam[j][1]==' ') atom.type=j+1;
    atoms.push_back( atom );
  }

  for( i=0; i<numbonds; i++) {
    fgets(s,80,f);
    sscanf(s,"%d %d %d\n", &ii, &bond.a, &bond.b);
    bond.o=0.1;
    if ((strstr(s,"SINGLE")!=NULL)) bond.o=0.1;
    if ((strstr(s,"DOUBLE")!=NULL)) bond.o=0.2;
    bonds.push_back(bond);
  }
  fclose(f);

  cread=1;
  bread=1;

  return 0;
}

/*
int ReadMolecule::XYZ (char* fname, int strn) {
  char s[256];
  FILE *f;
  double fl=0;
  struct atom atom;
  int strcur=0;
  int numatoms = 0;

  f=fopen(fname,"r");  if (f==NULL) return -1;

  do {

    fgets(s,250,f);
    sscanf(s,"%u %lf\n",&numatoms,&fl);

    if ( (numatoms == 0) || (fl != 0) ) {                                       // XYZ file without header
      fclose(f); 
      f=fopen(fname,"r");
      numatoms = 256;
    } else fgets(s,250,f);                                                      // comment

    for(int i=0; i<numatoms; i++) {

      s[0]=' '; fgets(s+1,250,f);
      if ( feof(f) != 0 ) break;
      if(s[1]=='*') { i--; continue; }                                          // comment
      sscanf(s,"%s %lf %lf %lf", &s, &atom.x, &atom.y, &atom.z);
      atom.type=0;                                                              // default
      if (atoi(s)>0) atom.type=atoi(s); else
      for (int j=0;j<at_namc;j++)
        if (s[0]==at_nam[j][0])
          if (s[1]==at_nam[j][1]) {atom.type=j+1; break; } else
            if (at_nam[j][1]==' ') atom.type=j+1;

      if ( cread == 1 ) atoms[i] = atom;
      else atoms.push_back( atom );

    }

    cread = 1;
    strcur++; 
    if(strcur==strn) { 
      fclose(f); 
      return 0; 
    }

  } while (feof(f)==0);

  fclose(f);

  tstrn=strcur-1;

  return 0;

}
*/

int ReadMolecule::XYZ (char* fname, int strn) {

  string s,ss;
  double fl=0;
  struct atom atom;
  int strcur=0;
  int numatoms = 0;

  ifstream f(fname);
  if ( !f ) return -1;

  while ( !f.eof() ) {

    getline(f,s);
    istringstream sin(s);
    sin >> numatoms >> fl;

    if ( (numatoms == 0) || (fl != 0) ) {                                       // XYZ file without header
      f.close(); 
      f.open(fname);
      numatoms = 256;
    } else getline(f,s);                                                      // comment

    for(int i=0; i<numatoms; i++) {

      getline(f,s);
      if ( f.eof() ) break;
      if(s.at(1)=='*') { i--; continue; }                                          // comment
  	  istringstream sin(s);
	  sin >> ss >> atom.x >> atom.y >> atom.z;
      
      atom.type=0;                                                              // default

  	  istringstream ssin(ss);
      ssin >> atom.type;
      if (atom.type == 0) {
        for (int j=0;j<at_namc;j++)
           if (ss.at(0)==at_nam[j][0])
             if (ss.size()==1) {atom.type=j+1; break; } else
               if (ss.at(1)==at_nam[j][1]) {atom.type=j+1; break; } else
                 if (at_nam[j][1]==' ') atom.type=j+1;
      }    

      if ( cread == 1 ) atoms[i] = atom;
      else atoms.push_back( atom );

    }

    cread = 1;
    strcur++; 
    if(strcur==strn) { 
      f.close(); 
      return 0; 
    }

  };

  f.close();

  tstrn=strcur-1;

  return 0;

}


void ReadMolecule::CenterMolecule() {                    // Centering molecule to (0,0,0)
  int i;

  float mx[2],my[2],mz[2];

        mx[0]=mx[1]=atoms[0].x; 
        my[0]=my[1]=atoms[0].y;
        mz[0]=mz[1]=atoms[0].z;

  for( i=1;i<atoms.size();i++ ) {
    if(atoms[i].x>mx[0]) mx[0]=atoms[i].x; else if(atoms[i].x<mx[1]) mx[1]=atoms[i].x;
    if(atoms[i].y>my[0]) my[0]=atoms[i].y; else if(atoms[i].y<my[1]) my[1]=atoms[i].y;
    if(atoms[i].z>mz[0]) mz[0]=atoms[i].z; else if(atoms[i].z<mz[1]) mz[1]=atoms[i].z;
  }

  mx[0]+=(mx[1]-mx[0])/2; 
  my[0]+=(my[1]-my[0])/2; 
  mz[0]+=(mz[1]-mz[0])/2;

  for( i=0;i<atoms.size();i++ ) { 
    atoms[i].x-=mx[0]; 
    atoms[i].y-=my[0]; 
    atoms[i].z-=mz[0]; 
  }
}

// Standart orientation of the Molecule
void ReadMolecule::ZOrientMolecule(int a0, int a1, int a2) {

  RowVector C(3);
  RowVector Z(3);

  int i;

  // Translate, atom 0 is 0,0,0
  Z(1) = atoms[a0].x;    Z(2) = atoms[a0].y;    Z(3) = atoms[a0].z;
  for( i=0; i<atoms.size(); i++ ) {
    C(1) = atoms[i].x;    C(2) = atoms[i].y;    C(3) = atoms[i].z;
    C -= Z;
    atoms[i].x = C(1);    atoms[i].y = C(2);    atoms[i].z = C(3);
  }

  SquareMatrix R(3),R2(3);

  // Rotate around z axis to zero x coordinate atom 1
  double leng = len( 0,0,0, atoms[a1].x,atoms[a1].y,0 );
  double sinz =  atoms[a1].x/leng;
  double cosz = -atoms[a1].y/leng;

  double newy1 = -sinz*atoms[a1].x+cosz*atoms[a1].y;
  double newx2 =  cosz*atoms[a2].x+sinz*atoms[a2].y;
  double newy2 = -sinz*atoms[a2].x+cosz*atoms[a2].y;

  R(1,1) = cosz;   R(1,2) = -sinz;   R(1,3) = 0.0;                              // z axis
  R(2,1) = sinz;   R(2,2) =  cosz;   R(2,3) = 0.0;
  R(3,1) =  0.0;   R(3,2) =   0.0;   R(3,3) = 1.0;

  // Rotate around x axis to zero y coordinate atom 1
  leng = len( 0,0,0, 0,newy1,atoms[a1].z );
  sinz = -newy1/leng;
  cosz =  atoms[a1].z/leng;

  newy2 = cosz*newy2+sinz*atoms[a2].z;

  R2(1,1) =  1.0;   R2(1,2) =   0.0;   R2(1,3) =   0.0;                         // x axis
  R2(2,1) =  0.0;   R2(2,2) =  cosz;   R2(2,3) = -sinz;
  R2(3,1) =  0.0;   R2(3,2) =  sinz;   R2(3,3) =  cosz;

  R *= R2;
  
  // Rotate around z axis to zero y coordinate atom 2
  leng = len( 0,0,0, newx2,newy2,0 );
  sinz =  newx2/leng;
  cosz = -newy2/leng;

  R2(1,1) = cosz;   R2(1,2) = -sinz;   R2(1,3) = 0.0;                           // z axis
  R2(2,1) = sinz;   R2(2,2) =  cosz;   R2(2,3) = 0.0;
  R2(3,1) =  0.0;   R2(3,2) =   0.0;   R2(3,3) = 1.0;

  R *= R2;
  
  for( i=0; i<atoms.size(); i++ ) {
    C(1) = atoms[i].x;    C(2) = atoms[i].y;    C(3) = atoms[i].z;
    C = C*R;
    atoms[i].x = C(1);    atoms[i].y = C(2);    atoms[i].z = C(3);
  }

}

// Determine the principal axes of inertia
void ReadMolecule::CalcPricipalAxes() {

  if (!cread) return;

  vector <struct atom>::iterator atom;
  int i;

  // find the center
  double sm=0,sx=0,sy=0,sz=0;
  for ( atom = atoms.begin(); atom!=atoms.end(); atom++ ) {
    sm+=(*atom).type;
    sx+=(*atom).x*(*atom).type;
    sy+=(*atom).y*(*atom).type;
    sz+=(*atom).z*(*atom).type;
  }

  // shift the origin
  double cx=sx/sm,cy=sy/sm,cz=sz/sm;
  for ( atom = atoms.begin(); atom!=atoms.end(); atom++ ) {
    (*atom).x-=cx;
    (*atom).y-=cy;
    (*atom).z-=cz;
  }

  // calculate the moments and products of inertia, ixx,iyy,izz and ixy,ixz,iyz 
  double ixx=0,iyy=0,izz=0,ixy=0,ixz=0,iyz=0;
  for( i=0; i<atoms.size(); i++ ) {
    ixx += atoms[i].type * ( sqr(atoms[i].y) + sqr(atoms[i].z) );
    iyy += atoms[i].type * ( sqr(atoms[i].x) + sqr(atoms[i].z) );
    izz += atoms[i].type * ( sqr(atoms[i].x) + sqr(atoms[i].y) );
    ixy -= atoms[i].type * atoms[i].x * atoms[i].y;
    iyz -= atoms[i].type * atoms[i].y * atoms[i].z;
    ixz -= atoms[i].type * atoms[i].x * atoms[i].z;
  }

  // inertia tensor I
  SymmetricMatrix I;
  I.resize(3);

  I(1,1) = ixx;  I(2,1) = ixy;  I(3,1) = ixz;
  I(1,2) = ixy;  I(2,2) = iyy;  I(3,2) = iyz;
  I(1,3) = ixz;  I(2,3) = iyz;  I(3,3) = izz;

  // diagonalize I  
  DiagonalMatrix I2;            // the principal moments
  SquareMatrix X;               // made up
  EigenValues(I, I2, X);

  for( i=0;i<atoms.size();i++) {
    double x=atoms[i].x, y=atoms[i].y, z=atoms[i].z;
    atoms[i].x = x*X(1,1) + y*X(2,1) + z*X(3,1);
    atoms[i].y = x*X(1,2) + y*X(2,2) + z*X(3,2);
    atoms[i].z = x*X(1,3) + y*X(2,3) + z*X(3,3);
  }

  I.release();
  I2.release();
  X.release();

}

