
#include <math.h>
#include <float.h>
#include <vector>

#include "geomath.h"

#define sqr(a)  ((a)*(a))

//#define atoms ReadMol->atoms

//extern ReadMoleculeBond *ReadMol;


double  getBondLen(atom a, atom b) {
   double Xa= a.x-b.x;
   double Ya= a.y-b.y;
   double Za= a.z-b.z;
   double q = Xa*Xa+Ya*Ya+Za*Za;
   return sqrt(q);
}


double  getAngle(atom a, atom c, atom b) {
  double 
   Xa= a.x-c.x,
   Ya= a.y-c.y,
   Za= a.z-c.z,
   Xb= b.x-c.x,
   Yb= b.y-c.y,
   Zb= b.z-c.z;
  double 
   r = Xa*Xb+Ya*Yb+Za*Zb,
   q = Xa*Xa+Ya*Ya+Za*Za,
   p = Xb*Xb+Yb*Yb+Zb*Zb;
  if (q <=0) {
//      printf("Atoms on single line: %d and %d\n",a+1,c+1);
      return 0;
  } else
  if (p <=0) {
//      printf("Atoms on single line: %d and %d\n",b+1,c+1);
      return 0;
  } else
  return 180.0/PI*acos(r/(sqrt(p)*sqrt(q)));
}


double  getDA(atom i1, atom i2, atom i3, atom i4) {  // from CCL
  double x1i,x2i,y1i,y2i,z1i,z2i,ux1,uy1,uz1,ux2,uy2,uz2,u1,u2,u,a,dihedr;

  x1i=i2.x-i1.x;
  y1i=i2.y-i1.y;
  z1i=i2.z-i1.z;
  x2i=i3.x-i2.x;
  y2i=i3.y-i2.y;
  z2i=i3.z-i2.z;

  ux1=y1i*z2i-z1i*y2i;
  uy1=z1i*x2i-x1i*z2i;
  uz1=x1i*y2i-y1i*x2i;

  x1i=i4.x-i3.x;
  y1i=i4.y-i3.y;
  z1i=i4.z-i3.z;
  
  ux2=z1i*y2i-y1i*z2i;
  uy2=x1i*z2i-z1i*x2i;
  uz2=y1i*x2i-x1i*y2i;

  u1=ux1*ux1+uy1*uy1+uz1*uz1;
  u2=ux2*ux2+uy2*uy2+uz2*uz2;
  u=u1*u2;

  if (u!=0.0) {
    a=(ux1*ux2+uy1*uy2+uz1*uz2)/sqrt(u);
    a=max(a,-1.0);
    a=min(a,1.0);
    dihedr=acos(a);
    if (ux1*(uy2*z2i-uz2*y2i)+
        uy1*(uz2*x2i-ux2*z2i)+
        uz1*(ux2*y2i-uy2*x2i) < 0.0) {dihedr =-dihedr;}
    return 180.0/PI*dihedr;
  } else {
//    printf("Error in dihedral (%d,%d,%d,%d)\n",i1+1,i2+1,i3+1,i4+1);
    return 0;
  }
}


void rotxyz (int a,int b,int c, double theta, float A[4][4]) {

  float    ct,st,B[9];

  ct = cos(PI*theta/180.0);
  st = sin(PI*theta/180.0);

  B[0] = A[0][0]*(a+ct*(1-a)) + A[0][1]*(  st*c)     + A[0][2]*( -st*b);
  B[1] = A[0][0]*( -st*c)     + A[0][1]*(b+ct*(1-b)) + A[0][2]*(  st*a);
  B[2] = A[0][0]*(  st*b)     + A[0][1]*( -st*a)     + A[0][2]*(c+ct*(1-c));

  B[3] = A[1][0]*(a+ct*(1-a)) + A[1][1]*(  st*c)     + A[1][2]*( -st*b);
  B[4] = A[1][0]*( -st*c)     + A[1][1]*(b+ct*(1-b)) + A[1][2]*(  st*a);
  B[5] = A[1][0]*(  st*b)     + A[1][1]*( -st*a)     + A[1][2]*(c+ct*(1-c));

  B[6] = A[2][0]*(a+ct*(1-a)) + A[2][1]*(  st*c)     + A[2][2]*( -st*b);
  B[7] = A[2][0]*( -st*c)     + A[2][1]*(b+ct*(1-b)) + A[2][2]*(  st*a);
  B[8] = A[2][0]*(  st*b)     + A[2][1]*( -st*a)     + A[2][2]*(c+ct*(1-c));

  A[0][0] = B[0];
  A[0][1] = B[1];
  A[0][2] = B[2];

  A[1][0] = B[3];
  A[1][1] = B[4];
  A[1][2] = B[5];

  A[2][0] = B[6];
  A[2][1] = B[7];
  A[2][2] = B[8];

}


void revmat(float A[4][4], float B[4][4]) {

  float
  det = A[0][0]*A[1][1]*A[2][2]
       -A[0][0]*A[1][2]*A[2][1]
       -A[1][0]*A[0][1]*A[2][2]
       +A[1][0]*A[0][2]*A[2][1]
       +A[2][0]*A[0][1]*A[1][2]
       -A[2][0]*A[0][2]*A[1][1];

  if (det != 0) {
    det = -1/det;
    B[0][0] = (-A[1][1]*A[2][2]+A[1][2]*A[2][1])*det;
    B[0][1] = ( A[0][1]*A[2][2]-A[0][2]*A[2][1])*det;
    B[0][2] = (-A[0][1]*A[1][2]+A[0][2]*A[1][1])*det;
    B[1][0] = ( A[1][0]*A[2][2]-A[1][2]*A[2][0])*det;
    B[1][1] = (-A[0][0]*A[2][2]+A[0][2]*A[2][0])*det;
    B[1][2] = ( A[0][0]*A[1][2]-A[0][2]*A[1][0])*det;
    B[2][0] = (-A[1][0]*A[2][1]+A[1][1]*A[2][0])*det;
    B[2][1] = ( A[0][0]*A[2][1]-A[0][1]*A[2][0])*det;
    B[2][2] = (-A[0][0]*A[1][1]+A[0][1]*A[1][0])*det;
  } else {
    B[0][0] = 1;  B[0][1] = 0;  B[0][2] = 0;
    B[1][0] = 0;  B[1][1] = 1;  B[1][2] = 0;
    B[2][0] = 0;  B[2][1] = 0;  B[2][2] = 1;
  }

}


// internal coordinates to cartesian
int int2car(atom i1, atom i2, atom i3, double R, double W, double T, double *X, double *Y, double *Z) {    // from Babel

  double cosph,sinph,costh,sinth,coskh,sinkh;
  double cosa,sina,cosd,sind;
  double xpd,ypd,zpd,xqd,yqd,zqd;
  double xa,ya,za;
  double rbc,xyb,yza,temp;
  double xpa,ypa,zqa;
  double xd,yd,zd;
  double x[5],y[5],z[5];

  x[1]=i1.x;  y[1]=i1.y;  z[1]=i1.z;
  x[2]=i2.x;  y[2]=i2.y;  z[2]=i2.z;
  x[3]=i3.x;  y[3]=i3.y;  z[3]=i3.z;

  double
    dist  = R,
    angle = W * DEG_TO_RAD,
    dihed = T * DEG_TO_RAD;

  const int i=4, na=3, nb=2, nc=1;

  double
    xb = x[nb] - x[na],
    yb = y[nb] - y[na],
    zb = z[nb] - z[na];
    
  rbc = xb*xb + yb*yb + zb*zb;
  if ( rbc < 0.0001 ) return 1;
  rbc = 1.0/sqrt(rbc);
    
  cosa = cos(angle);
  sina = sin(angle);
    
  if( fabs(cosa) >= 0.999999 ) {       // Colinear
    temp = dist*rbc*cosa;
    x[i] = x[na] + temp*xb;
    y[i] = y[na] + temp*yb;
    z[i] = z[na] + temp*zb;
  } else {
    xa = x[nc] - x[na];
    ya = y[nc] - y[na];
    za = z[nc] - z[na];
    
    sind = -sin(dihed);
    cosd =  cos(dihed);
    
    xd = dist*cosa;
    yd = dist*sina*cosd;
    zd = dist*sina*sind;
    
    xyb = sqrt(xb*xb + yb*yb);

    int flag = 0;
    
    if( xyb < 0.1 ) {                 // Rotate about y-axis
        temp = za; za = -xa; xa = temp;
        temp = zb; zb = -xb; xb = temp;
        xyb = sqrt(xb*xb + yb*yb);
        flag = 1;
    }   
    
    costh = xb/xyb;
    sinth = yb/xyb;
    xpa = costh*xa + sinth*ya;
    ypa = costh*ya - sinth*xa;
    
    sinph = zb*rbc;
    cosph = sqrt(1.0 - sqr(sinph));
    zqa = cosph*za  - sinph*xpa;
    
    yza = sqrt(ypa*ypa + zqa*zqa);
    
    if( yza > 1.0E-10 ) {   
        coskh = ypa/yza;
        sinkh = zqa/yza;
    } else { 
        coskh = 1.0;
        sinkh = 0.0;
    }
    
    ypd = coskh*yd  - sinkh*zd;
    zpd = coskh*zd  + sinkh*yd;
    xpd = cosph*xd  - sinph*zpd;
    zqd = cosph*zpd + sinph*xd;
    xqd = costh*xpd - sinth*ypd;
    yqd = costh*ypd + sinth*xpd;
    
    if( flag ) {                                        // Rotate about y-axis
        *X = x[na] - zqd;
        *Y = y[na] + yqd;
        *Z = z[na] + xqd;
    } else {  
        *X = x[na] + xqd;
        *Y = y[na] + yqd;
        *Z = z[na] + zqd;
    }
  }

  return 0;
}


// Translate the molecule relative to the screen
void shiftxyz(atomvector atoms, double sx, double sy, double sz, float roma[4][4], float romb[4][4]) {

  revmat(roma,romb);                                    // calculate A^(-1), where A - rotation matrix

  double
     tx = romb[0][0]*sx + romb[1][0]*sy + romb[2][0]*sz,
     ty = romb[0][1]*sx + romb[1][1]*sy + romb[2][1]*sz,
     tz = romb[0][2]*sx + romb[1][2]*sy + romb[2][2]*sz;

  for ( int i=0; i<atoms.size(); i++) {
     atoms[i].x += tx;
     atoms[i].y += ty;
     atoms[i].z += tz;
  }

}


void cpymat(float A[4][4], float B[4][4]) {
  for ( int i=0; i<3; i++) {
     for ( int j=0; j<3; j++) {
        B[i][j] = A[i][j];
     }
  }
}


// rotate the vector
void rotvec(float V[3], float A[4][4]) {

  double sx = V[0], sy = V[1], sz = V[2];
  V[0] = A[0][0]*sx + A[1][0]*sy + A[2][0]*sz,
  V[1] = A[0][1]*sx + A[1][1]*sy + A[2][1]*sz,
  V[2] = A[0][2]*sx + A[1][2]*sy + A[2][2]*sz;

}


// rotate the whole system
void rotall(atomvector atoms, float A[4][4]) {

  for ( int i=0; i<atoms.size(); i++) {
     double sx = atoms[i].x, sy = atoms[i].y, sz = atoms[i].z;
     atoms[i].x = A[0][0]*sx + A[1][0]*sy + A[2][0]*sz,
     atoms[i].y = A[0][1]*sx + A[1][1]*sy + A[2][1]*sz,
     atoms[i].z = A[0][2]*sx + A[1][2]*sy + A[2][2]*sz;
  }

}


// distance between point (x3,y3,z3) and line (x1,y1,z1)-(x2,y2,z2)
double LinePointDistance(double x1, double y1, double z1,
                         double x2, double y2, double z2,
                         double x3, double y3, double z3  ) {
   double x4 = (-2*x3*x2*x1+x2*x2*x3+x1*x1*x3-y1*x1*y2-z1*x1*z2+y2*y2*x1-
                 y1*x2*y3-y2*y1*x2+y2*y3*x2-z1*x2*z3+z2*z3*x2-z2*z1*x2-z3*z2*x1+
                 x1*z3*z1-x1*y3*y2+x1*y3*y1+z1*z1*x2+z2*z2*x1+y1*y1*x2)/
                (y1*y1+z2*z2+z1*z1-2*x2*x1-2*y2*y1-2*z2*z1+y2*y2+x2*x2+x1*x1);
   double y4 = (-x4 * y2 + x4 * y1 + x1 * y2 - y1 * x2) / (-x2 + x1);
   double z4 = (-x4 * z2 + x4 * z1 + x1 * z2 - z1 * x2) / (-x2 + x1);
   double r2=sqr(x4-x3)+sqr(y4-y3)+sqr(z4-z3);
   return sqrt(r2);
}


// distance between point (x3,y3,z3) and line segment (x1,y1,z1)-(x2,y2,z2)
double SegmentPointDistance(double x1, double y1, double z1,
                            double x2, double y2, double z2,
                            double x3, double y3, double z3  ) {
   double x4 = (-2*x3*x2*x1+x2*x2*x3+x1*x1*x3-y1*x1*y2-z1*x1*z2+y2*y2*x1-
                 y1*x2*y3-y2*y1*x2+y2*y3*x2-z1*x2*z3+z2*z3*x2-z2*z1*x2-z3*z2*x1+
                 x1*z3*z1-x1*y3*y2+x1*y3*y1+z1*z1*x2+z2*z2*x1+y1*y1*x2)/
                (y1*y1+z2*z2+z1*z1-2*x2*x1-2*y2*y1-2*z2*z1+y2*y2+x2*x2+x1*x1);
   double y4 = (-x4 * y2 + x4 * y1 + x1 * y2 - y1 * x2) / (-x2 + x1);
   double z4 = (-x4 * z2 + x4 * z1 + x1 * z2 - z1 * x2) / (-x2 + x1);   // x4,y4,z4 - point on the line (x1,y1,z1)-(x2,y2,z2)
                                                                        // intersection with perpendicular line from (x3,y3,z3)
   double r1 = sqrt(sqr(x4-x1)+sqr(y4-y1)+sqr(z4-z1));                  // point to segment begin distance 
   double r2 = sqrt(sqr(x4-x2)+sqr(y4-y2)+sqr(z4-z2));                  // point to segment end   distance 
   double r3 = sqrt(sqr(x2-x1)+sqr(y2-y1)+sqr(z2-z1));                  // segment length

   double rr;                                                           // distance from point to segment
   
   if ((r1+r2-r3)<0.00001) {                                            // point is inside the segment
      rr = sqrt(sqr(x4-x3)+sqr(y4-y3)+sqr(z4-z3));                      // calculate distance to line
   } else 
      if (r1<r2) {                                                      // point closer to begin of segment
        rr = sqrt(sqr(x3-x1)+sqr(y3-y1)+sqr(z3-z1));                    // distance to begin
      } else {                                                          // point closer to end   of segment
        rr = sqrt(sqr(x3-x2)+sqr(y3-y2)+sqr(z3-z2));                    // distance to end
      }
   
   return rr;
}

