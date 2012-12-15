
#include <cmath>

#define PI 3.14159265358979323846

#ifdef WIN32
#include <windows.h>
#include <SHLWAPI.H>
#endif

#include <GL/gl.h>
#include <GL/glu.h>

#define X 0.525731112119133606
#define Z 0.850650808352039932

#define SphereList    1
#define CylinderList  2

static float vdata[12][3] = {
  {-X,0.0,Z}, {X,0.0,Z}, {-X,0.0,-Z}, {X,0.0,-Z},
  {0.0,Z,X}, {0.0,Z,-X}, {0.0,-Z,X}, {0.0,-Z,-X},
  {Z,X,0.0}, {-Z,X,0.0}, {Z,-X,0.0}, {-Z,-X,0.0}
};

static int tindices[20][3] = {
  {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
  {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
  {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
  {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
};

// normalize a vector of non-zero length
void normalize(float v[3]) {
  float d = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] /= d; v[1] /= d; v[2] /= d;
}

/* recursively subdivide face `depth' times */
/* and draw the resulting triangles */
void subdivide(GLfloat v1[3], GLfloat v2[3], GLfloat v3[3], int depth)
{
  GLfloat v12[3], v23[3], v31[3];
  int i;

  if (depth == 0) {
    glBegin(GL_TRIANGLES);
      glNormal3fv(v1);   glVertex3fv(v1);      
      glNormal3fv(v2);   glVertex3fv(v2);      
      glNormal3fv(v3);   glVertex3fv(v3);      
    glEnd(); 
    return;
  }

  /* calculate midpoints of each side */
  for (i = 0; i < 3; i++) {
    v12[i] = (v1[i]+v2[i])/2.0;
    v23[i] = (v2[i]+v3[i])/2.0;
    v31[i] = (v3[i]+v1[i])/2.0;
  }

  /* extrude midpoints to lie on unit sphere */
  normalize(v12);
  normalize(v23);
  normalize(v31);

  /* recursively subdivide new triangles */
  subdivide(v1 , v12, v31, depth-1);
  subdivide(v2 , v23, v12, depth-1);
  subdivide(v3 , v31, v23, depth-1);
  subdivide(v12, v23, v31, depth-1);
}

int olddepth = -1;

void CreateSphereList(int depth) {
  glNewList(SphereList, GL_COMPILE);
  for (int i = 0; i < 20; i++) {
    subdivide(&vdata[tindices[i][0]][0],
              &vdata[tindices[i][1]][0],
	          &vdata[tindices[i][2]][0],
	          depth);
  }
  glEndList();
}

void Sphere(GLdouble radius, GLint slices, GLint stacks ) {
  int depth = slices/16;

  if (depth != olddepth) {
     CreateSphereList(depth);
  }

  glScalef(radius,radius,radius);
  glCallList(SphereList);
}


/*
#define MAXsl 140

GLdouble sinez[MAXsl], cosez[MAXsl];

int SphereLastSlices = 0;

void Sphere(GLdouble radius, GLint slices, GLint stacks ) {
   GLdouble x, y, z;
   GLint i, j, imax;

   if(slices!=SphereLastSlices) {
     for (j=0;j<=slices;j++) {
       sinez[j]=sin((GLdouble)j * (PI+PI) / (GLdouble) slices);
       cosez[j]=cos((GLdouble)j * (PI+PI) / (GLdouble) slices);
     }
     SphereLastSlices=slices;
   }

   glBegin( GL_TRIANGLE_FAN );
   glNormal3d( 0.0, 0.0, 1.0 );
   glVertex3d( 0.0, 0.0, radius );
   for (j=0;j<=slices;j++) {
      x = -sinez[j] * sinez[1];
      y =  cosez[j] * sinez[1];
      z =  cosez[1];
      glNormal3d( x, y, z );
      glVertex3d( x*radius, y*radius, z*radius );
   }
   glEnd();

   // draw intermediate stacks as quad strips
   for (i=1,imax = stacks/2;i<imax;i++) {
      glBegin( GL_QUAD_STRIP );
      for (j=0;j<=slices;j++) {
         x = -sinez[j] * sinez[i];
         y =  cosez[j] * sinez[i];
         z =  cosez[i];
         glNormal3d( x, y, z );
         glVertex3d( x*radius, y*radius, z*radius );
         x = -sinez[j] * sinez[i+1];
         y =  cosez[j] * sinez[i+1];
         z =  cosez[i+1];
         glNormal3d( x, y, z );
         glVertex3d( x*radius, y*radius, z*radius );
      }
      glEnd();
   }

   glBegin( GL_TRIANGLE_FAN );
   glNormal3d( 0.0, 0.0, -1.0 );
   glVertex3d( 0.0, 0.0, -radius );
   
   for (i=slices/2-1, j=slices;j>=0;j--) {
      x = -sinez[j] * sinez[i];
      y =  cosez[j] * sinez[i];
      z =  cosez[i];
      glNormal3d( x, y, z );
      glVertex3d( x*radius, y*radius, z*radius );
   }
   glEnd();
}*/

int oldslices = -1;

void Cylinder(double baseRadius, double topRadius, 
              double height,            int slices) {

  if (slices < 3) slices = 3;

  if (slices != oldslices) {                // create new cylinder list

       GLdouble x1c, y1c, x2c, y2c;

       double da = 2.0*PI / slices;
       double dr = (topRadius-baseRadius);
       double nz = -dr;
   
       glNewList(CylinderList, GL_COMPILE);

       for (int i=0; i<slices; i++) {
         x1c = -sin( i   *da);
         y1c =  cos( i   *da);
         x2c = -sin((i+1)*da);
         y2c =  cos((i+1)*da);

         double r = baseRadius;
         glBegin( GL_QUAD_STRIP );
         glNormal3d( x1c  , nz ,   y1c);
         glVertex3d( x1c*r, 0  , y1c*r);
         glNormal3d( x2c  , nz ,   y2c);
         glVertex3d( x2c*r, 0  , y2c*r);
         r += dr;
         glNormal3d( x1c  , nz ,   y1c);
         glVertex3d( x1c*r, 1.0, y1c*r);
         glNormal3d( x2c  , nz ,   y2c);
         glVertex3d( x2c*r, 1.0, y2c*r);
         glEnd();
       }
      
       glEndList();
  }

  glScalef(1.0, height, 1.0);
  glCallList(CylinderList);

}
