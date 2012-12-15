
using namespace std;

#include "cr_bond.h"

const double PI=3.14159265358979323846;
const double DEG_TO_RAD = PI/180.0;

double  getBondLen(atom a,  atom b);
double  getAngle  (atom a,  atom c,  atom b);
double  getDA     (atom i1, atom i2, atom i3, atom i4);

int     int2car   (atom i1, atom i2, atom i3, 
                   double R, double W, double T, 
                   double *X, double *Y, double *Z);

void    rotxyz    (int a,int b,int c, double theta, float A[4][4]);
void    revmat    (float A[4][4], float B[4][4]);                       // inverse A to B
void    cpymat    (float A[4][4], float B[4][4]);                       // copy A to B
void    shiftxyz  (atomvector atoms, double, double, double, float, float);
void    rotvec    (float V[3], float A[4][4]);
void    rotall    (atomvector atoms, float A[4][4]);

double LinePointDistance(double x1, double y1, double z1,
                         double x2, double y2, double z2,
                         double x3, double y3, double z3  );

double SegmentPointDistance(double x1, double y1, double z1,
                            double x2, double y2, double z2,
                            double x3, double y3, double z3  );
