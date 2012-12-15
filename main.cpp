#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <ctime>

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Hor_Slider.H>
#include <FL/Fl_Roller.H>
#include <FL/Fl_Toggle_Button.H>
#include <FL/Fl_Menu_Button.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Help_Dialog.H>
#include <FL/Fl_Progress.H>
#include <FL/Fl_Multiline_Output.H>
#include <FL/fl_draw.H>
#include <FL/gl.h>
#include <GL/glu.h>

#include "geomath.h"
#include "gl_param.h"
#include "atom_rads.h"
#include "graphprim.h"
#include "lbfgs/lbfgs.h"

#define atm(a,x) atoms[bonds[i].a-1].x
#define ax1 atm(a,x)
#define ax2 atm(b,x)
#define ay1 atm(a,y)
#define ay2 atm(b,y)
#define az1 atm(a,z)
#define az2 atm(b,z)

#define len(x,y,z,xx,yy,zz) sqrt(sqr((x)-(xx))+sqr((y)-(yy))+sqr((z)-(zz)))
#define sqr(a)  ((a)*(a))

const int 
            tw = 400,                 // message window width
            wh = 600, ww = wh + tw,   // main window height and width
            sw = 20,                  // width or height of axes rotation sliders
            sr = 30,                  // height of zoom roller
            ds = 5,                   // distance between elements
            rw = 60;                  // roller width

const char * AboutText = 
            "<HTML>\n"
            "<HEAD>\n"
            "<TITLE>ViewMol3D</TITLE>\n"
            "</HEAD>\n"
            "<BODY BGCOLOR='#ffffff'>\n"
            "<h1 align=center>ViewMol3D</h1>\n"
            "<h2 align=center>A 3D OpenGL viewer for molecular structures<br>\n"
            "from the output of quantum chemistry calculations.<br></h2>\n"
            "<p align=left><font size=4>ViewMol3D uses OpenGL (OpenGL is a trademark of Silicon Graphics," 
            "Inc.) 3D graphical system to render molecules as wire frame, sticks, ball-and-sticks and CPK models. The ViewMol3D can draw" 
            "molecules models from output of several quantum chemistry programs.</font></p>"
            "<p align=left><font color=blue size=5>Possibilities</font>"
            "<font size=4>"
            "<ul>"
            "  <li><p align=left>Showing the geometry of a molecule</p></li>"
            "  <li><p align=left>Tracing a geometry optimization or a MD trajectory</p></li>"
            "  <li><p align=left>Showing normal vibrations of a molecule as arrows</p></li>"
            "  <li><p align=left>Showing forces acting on each atom in a selected configuration</p></li>"
            "  <li><p align=left>Saving all generated pictures as BMP/PNG file.</p></li>"
            "</ul>" 
            "</font>"
            "<p align=left><font color=blue size=5>File Formats Supported</font>"
            "<font size=4>"
            "<ul>"
            "  <li><p align=left>Gaussian 9x/03 OUT files</p></li>"
            "  <li><p align=left>GAMESS(US)/PC GAMESS OUT files</p></li>"
            "  <li><p align=left>MOPAC/AMPAC OUT files</p></li>"
            "  <li><p align=left>XYZ files with Cartesian coordinates</p></li>"
            "  <li><p align=left>Tripos Alchemy MOL files</p></li>"
            "</ul>"
            "</font>"
//            "<H2>ViewMol3D</H2>\n"
//            "<P>is a program for visualization of molecular models.</P>\n"
            "</BODY>\n"
            "</HTML>";

const int                                         // selection levels
            at_level = 1000,                      // atoms
            bn_level = 3000,                      // bonds
            gr_level = 4000;                      // gradients

mindo3 *mcalc;
atomvector_ehm atoms_ehm;


class MyGlWindow : 
	public Fl_Gl_Window, 			// OpenGL window and interface elements
	public ReadMoleculeBond, 		// Molecule specification and reading from files
	public LBFGS				    // L-BFGS optimization
{

  public:

    // Name of the input file
    char fname[256];

    // global variables
    int
       trad,                   // true if van-der-vaals radii
       bond,                   // is there bonds or not
       persp,                  // perpective
       LowQual,                // bonds are wires
       dline,                  // bonds are lines
       chg,                    // show atomic charges
       present,                // presentation mode
       monum,                  // MO number
       atnum,                  // show atomic numbers
       viewir,                 // show freq vecs
       frcur,                  // cursor position for printing
       viewgr,                 // show grad vecs
       strnum,                 // number of structure to show
       xaccl,                  // graphic accelerator is used
       screenNum,              // number of saved screens
       irc,                    // show IRC graph
       clrbnd,                 // flag of colored bonds
       pbuff,                  // pixel buffers are absent
       mrkatm,                 // show marked atoms
       DrawMO,                 // viaualize molecular orbitals
       DrawAxes,               // xyz axes
       DrawSurf                // surface
    ;

    int ap1, ap2, ap3, ap4;    // marked atoms numbers

    int x_select,y_select;     // coordinates of the selection
    int select_mode;           // default mode GL_RENDER, not GL_SELECT
    int select_hit;            // which object was select 

    int slnormalsphere, slnormalcylinder, slpresentsphere, slpresentcylinder;
    double quality, old_quality;

    double old_angle_x, angle_x;
    double old_angle_y, angle_y;
    double old_angle_z, angle_z;

    int numorb;
    double Energy, RMS;
    int funcalls;
    
    Fl_Window       *win;
    Fl_Menu_Bar     *menubar;
    Fl_Text_Buffer  *buff;
    Fl_Text_Display *disp;
    Fl_Roller       *zoom;
    Fl_Hor_Slider   *slider_x;
    Fl_Slider       *slider_y;
    Fl_Slider       *slider_z;
    Fl_Button       *get_len;
    Fl_Button       *get_ang;
    Fl_Button       *get_dih;
    Fl_Button       *get_zmt;
    Fl_Help_Dialog  *about;

    int ReadInputFile() {
          
              int error=0;
              char  xstr[64];
          
              strcpy (xstr, &fname[strlen(fname)-3]);
              if ((strncmp (xstr, "out", 3) == 0) || (strncmp (xstr, "OUT", 3) == 0) ||
                  (strncmp (xstr, "gms", 3) == 0) || (strncmp (xstr, "GMS", 3) == 0) ||
                  (strncmp (xstr, "log", 3) == 0) || (strncmp (xstr, "LOG", 3) == 0) ||
                  (strncmp (xstr, "mop", 3) == 0) || (strncmp (xstr, "MOP", 3) == 0))   {
                      if (OutFile (fname,strnum) != 0) {
    //                        sprintf(lpInfo,"Can't open the GAMESS/G98/MOPAC OUT file '%s'", fname);
                              error=1;
                      }
              } else if ( (strncmp (xstr, "mol", 3) == 0) || (strncmp (xstr, "MOL", 3) == 0) ) {
                      if (ALCHEMY (fname) != 0) {
    //                        sprintf(lpInfo,"Can't open the ALCHEMY MOL file '%s'", fname);
                              error=1;
                      } else dline=0;
              } else if ( (strncmp (xstr, "xyz", 3) == 0) || (strncmp (xstr, "XYZ", 3) == 0) ) {
    //          } else if (strncasecmp (xstr, "xyz", 3) == 0) {
                      if (XYZ (fname,strnum) != 0) {
    //                        sprintf(lpInfo,"Can't open the XYZ file '%s'", fname);
                              error=1;
                      } else dline=0;
              } else {
    //                sprintf(lpInfo,"I don't know what kind of this file:'%s'. Please give me a suitable extention (MOL/OUT/LOG/GMS/MOP)",fname);
                      error=1;
              }
          
              if ( error != 0 ) {
    //                MessageBox(NULL, lpInfo, "Error", MB_ICONERROR | MB_OK);
                      return FALSE;
              } else return TRUE;
          
    }


    double maxcor, x0,y0,xm,ym,z0,zm;
    float  roma[4][4], romb[4][4];

    void CalcMaxCor() {

//       if ( ReadMol == NULL ) return;
       if ( atoms.size() <= 0 ) return;

       float xx=fabs(atoms[0].x);
           for(int i=0;i<atoms.size();i++) {
                 if(fabs(atoms[i].x)>xx) xx=fabs(atoms[i].x);
                 if(fabs(atoms[i].y)>xx) xx=fabs(atoms[i].y);
                 if(fabs(atoms[i].z)>xx) xx=fabs(atoms[i].z);
           }
           maxcor=xx+xx;

       x0=-maxcor/1.2; xm=maxcor/1.2;
       y0=-maxcor/1.2; ym=maxcor/1.2;
       z0=-maxcor/1.2; zm=maxcor/1.2;
    }

    
    void SetOrtho() {
        
        if (select_mode == GL_SELECT) {
            int vp[4];                                                      // view port
            glGetIntegerv(GL_VIEWPORT, vp);
            glViewport (0, 0, w(), h());
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluPickMatrix(x_select, h()-y_select, 1, 1, vp);
        } else {
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
        }
        
        double ww=w(); double hh=h();
        double wh  = ww/hh;
        double mwh = maxcor*wh;
        
        if (!persp) {
              if (ww <= hh) 
                 glOrtho (-maxcor, maxcor,  -mwh,    mwh,    -maxcor, maxcor);
              else 
                 glOrtho (-mwh,    mwh,     -maxcor, maxcor, -maxcor, maxcor);
        } else {
                gluPerspective (50.0, wh, 1.0, 50.0);
        }
        
        glMatrixMode (GL_MODELVIEW);
        glLoadIdentity ();
        if (persp) glTranslated (0, 0, -2 * maxcor);  // viewing transform
    }

    
    void InitGL() {

        glClearColor(0., 0., 0., 0.);
        glViewport (0, 0, w(), h());
//        if (persp) glTranslated (0, 0, -2 * maxcor);        // viewing transform

        SetOrtho();

        glEnable(GL_AUTO_NORMAL);

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glClearDepth(1.0);

        float
            light_position1[4]= {float(-maxcor/2),float(-maxcor/2),float(maxcor  ),0},
            light_position2[4]= {float( maxcor/2),float( maxcor/2),float(maxcor/2),0},
            light_position3[4]= {float( 0.0     ),float( 0.0     ),float(maxcor  ),0};

        const float
            lm_ambient1          [4] =  { 0.2f, 0.2f, 0.2f, 1.0f },
            light_color1         [4] =  { 1.0f, 0.8f, 0.7f, 1.0f },
            light_color2         [4] =  { 0.7f, 0.9f, 1.0f, 1.0f },
            light_color3         [4] =  { 0.2f, 0.2f, 0.2f, 1.0f },
            light_color_diffuse1 [4] =  { 0.6f, 0.6f, 0.6f, 1.0f },
            light_color_diffuse2 [4] =  { 0.6f, 0.6f, 0.6f, 1.0f },
            light_color_diffuse3 [4] =  { 0.6f, 0.6f, 0.6f, 1.0f };
/*        
        float v1[4],v2[4],v3[4],v4[4],v5[3],v6,v7,v8,v9,v10; int v11;
        glGetLightfv(GL_LIGHT0,GL_AMBIENT,v1);                      printf("1 %f %f %f %f\n",v1[0],v1[1],v1[2],v1[3]);
        glGetLightfv(GL_LIGHT0,GL_DIFFUSE,v2);                      printf("2 %f %f %f %f\n",v2[0],v2[1],v2[2],v2[3]);
        glGetLightfv(GL_LIGHT0,GL_SPECULAR,v3);                     printf("3 %f %f %f %f\n",v3[0],v3[1],v3[2],v3[3]);
        glGetLightfv(GL_LIGHT0,GL_POSITION,v4);                     printf("4 %f %f %f %f\n",v4[0],v4[1],v4[2],v4[3]);
        glGetLightfv(GL_LIGHT0,GL_SPOT_DIRECTION,v5);               printf("5 %f %f %f   \n",v5[0],v5[1],v5[2]);
        glGetLightfv(GL_LIGHT0,GL_SPOT_EXPONENT,&v6);               printf("6 %f         \n",v6);
        glGetLightfv(GL_LIGHT0,GL_SPOT_CUTOFF,&v7);                 printf("7 %f         \n",v7);
        glGetLightfv(GL_LIGHT0,GL_CONSTANT_ATTENUATION,&v8);        printf("8 %f         \n",v8);
        glGetLightfv(GL_LIGHT0,GL_LINEAR_ATTENUATION,&v9);          printf("9 %f         \n",v9);
        glGetLightfv(GL_LIGHT0,GL_QUADRATIC_ATTENUATION,&v10);      printf("A %f         \n",v10);
*/
        glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_color1);
//        glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_color_diffuse1);

        glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
        glLightfv(GL_LIGHT1, GL_SPECULAR, light_color2);
//        glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_color_diffuse2);

        glLightfv(GL_LIGHT2, GL_POSITION, light_position3);
        glLightfv(GL_LIGHT2, GL_SPECULAR, light_color3);
//        glLightfv(GL_LIGHT2, GL_DIFFUSE,  light_color_diffuse3);

        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lm_ambient1);

        glEnable(GL_LIGHT0);
        glEnable(GL_LIGHT1);
//        glEnable(GL_LIGHT2);

        glEnable(GL_LIGHTING);

        if (present) glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST); 
    
        glEnable(GL_NORMALIZE);

        glShadeModel(GL_SMOOTH);

        const float
            mat_ambient1    [4] =  { 1.0f, 1.0f, 1.0f, 1.0f },
            mat_diffuse1    [4] =  { 0.6f, 0.6f, 0.6f, 1.0f },
            mat_specular1   [4] =  { 0.8f, 0.8f, 0.8f, 0.8f },
            mat_emission1   [4] =  { 0.0f, 0.0f, 0.0f, 1.0f };

        glMaterialfv(GL_FRONT, GL_AMBIENT  , mat_ambient1 );
        glMaterialfv(GL_FRONT, GL_DIFFUSE  , mat_diffuse1 );
        glMaterialfv(GL_FRONT, GL_SPECULAR , mat_specular1);
        glMaterialfv(GL_FRONT, GL_EMISSION , mat_emission1);
        glMaterialf (GL_FRONT, GL_SHININESS,  40.f ); 

    }


    void drawmark (float color[4], float angle, float rad, int sl) {
         glMaterialfv (GL_FRONT, GL_DIFFUSE, color); 
         glRotatef(angle, 1.0, 1.0, 1.0);
//         glutWireSphere(rad,sl/2,sl/2);
         glTranslatef(0.0,-rad/10.0,0.0);                                    Cylinder(rad,rad,rad/5.0,sl);
         glTranslatef(0.0,rad/10.0,-rad/10.0); glRotatef(90, 1.0, 0.0, 0.0); Cylinder(rad,rad,rad/5.0,sl);
         glTranslatef(rad/10.0,rad/10.0,0.0);  glRotatef(90, 0.0, 0.0, 1.0); Cylinder(rad,rad,rad/5.0,sl);

    }

    
    void print_label(int n, int i, float c[3]) {
        char p[32];
        if (n>0)
          sprintf(p,"%u[%u]",i+1,n); 
        else
          sprintf(p,"%u",i+1); 
        rotvec (c, romb);
        glRasterPos3f( atoms[i].x+c[0],                    // label position
                       atoms[i].y+c[1],                    // plus atom radius
                       atoms[i].z+c[2]  );
        if (p[0]!=0) gl_draw(p, strlen(p));
    }
    
    // DRAW METHOD
    void draw() {

        // First time: init viewport, etc.
        if ( !valid() ) {
            valid(1);
            CalcMaxCor();
            InitGL();
            char msg[256];
            sprintf(msg,"OpenGL initialized:\n"
                            "     Vendor   : %s\n"
                            "     Renderer : %s\n"
                            "     Version  : %s\n"
                            " GLU Version  : %s\n",
                            glGetString  (GL_VENDOR  ),
                            glGetString  (GL_RENDERER),
                            glGetString  (GL_VERSION ),
                            gluGetString (GLU_VERSION));
            buff->text(msg);
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                         // Clear screen

        glPushMatrix();

        glMultMatrixf((float*)roma);                                                // rotate the system

        glPushMatrix();

        // draw atoms
        int sl = bond ?
              sl = present ? slpresentcylinder : slnormalcylinder             // No atoms, just bonds
        :
              sl = present ? slpresentsphere   : slnormalsphere
        ;

        int i=0;

        for ( vector <struct atom>::iterator it = atoms.begin(); it!=atoms.end(); it++, i++) {

              glPushMatrix ();

              float rad;

              glTranslated (atoms[i].x, atoms[i].y, atoms[i].z);

              if (select_mode == GL_SELECT) glLoadName(at_level+i);

              // set material
              if (chg) {                                                            // atoms charge
                      if(atoms[i].c>0) glMaterialfv (GL_FRONT, GL_DIFFUSE, PosChg);
                      else             glMaterialfv (GL_FRONT, GL_DIFFUSE, NegChg);
                      rad=fabs(atoms[i].c)*2;
              } else {
                      if (bond) {                                                   // no atoms 
                        if(clrbnd) glMaterialfv (GL_FRONT, GL_DIFFUSE, at_def[atoms[i].type-1]);
                        else           glMaterialfv (GL_FRONT, GL_DIFFUSE, bondmat);
                        rad=0;                                                      // if(dline) - bonds are linear
                        if(!dline)                                                  // if bonds are cylindrical
                              for ( vector <struct bond>::iterator it = bonds.begin(); it!=bonds.end(); it++) { // Radius of a biggest bonds cylinder
                                if ((((*it).a-1==i)||((*it).b-1==i))&&((*it).o>rad)) rad=(*it).o;
                              }    
                      } else {                                                      // Atoms are drawn
                        glMaterialfv (GL_FRONT, GL_DIFFUSE, at_def[atoms[i].type-1]);
                        if (trad==0) rad=at_rad [atoms[i].type-1]/2;                // Normal radii
                        else         rad=at_radv[atoms[i].type-1];                  // Van-der-Waals radii
                      }
              }

              Sphere(rad,sl,sl);                     // Draw the Atom

              if ( mrkatm ) {
                if ( (i==ap1) || (i==ap2) || (i==ap3) || (i==ap4) ) {
                       rad+=0.02f;
                       if ( i == ap1 ) drawmark(White, 0,rad,sl);
                       if ( i == ap2 ) drawmark(Red  ,30,rad,sl);
                       if ( i == ap3 ) drawmark(Green,60,rad,sl);
                       if ( i == ap4 ) drawmark(Blue ,90,rad,sl);
                }
              }

              glPopMatrix ();
        }

        // draw atomic numbers
        if (atnum) {

            revmat(roma,romb);                                              // find matrix inverse to the rotation matrix

            int h = int(double(at_rad[0])*double(w())/(maxcor*3));          // font should be apprx equal to the hydrogen atom radius
            if ( h<8 ) h=8;                                                 // but not too small
            gl_font(1, h);
            glMaterialfv (GL_FRONT, GL_DIFFUSE, Blue);                      // blue material for label

            for ( int i=0; i<atoms.size(); i++ ) {

                float rad;

                if (trad==0) rad=at_rad [atoms[i].type-1]/2;       // Normal radii
                else         rad=at_radv[atoms[i].type-1];         // Van-der-Waals radii

                float v[3] = {0.0f, 0.0f, 1.0f};                            // create normal
                rotvec (v, romb);                                           // perpendicular to the screen

                glNormal3fv(v);                                             // the normal is perpendicular in screen direction
                
                float c[3] = {0.0f, 0.0f, rad+0.1f};                        // coordinates of the label
                if ( mrkatm ) {
                    if ( i == ap1 ) { print_label(1,i,c); c[1]-=0.25; }
                    if ( i == ap2 ) { print_label(2,i,c); c[1]-=0.25; } 
                    if ( i == ap3 ) { print_label(3,i,c); c[1]-=0.25; }
                    if ( i == ap4 ) { print_label(4,i,c); c[1]-=0.25; }
                    if ( (i!=ap1) && (i!=ap2) && (i!=ap3) && (i!=ap4) ) {              
                                      print_label(0,i,c);             
                    }
                } else {              
                                      print_label(0,i,c);             
                }

           }
        }

        // Draw Bonds
        if (trad==0) {                                                      // if not CPK model, which doesn't need bonds

            if(dline) {                                                     // wire bonds

                float v1[4],v2[4],v3[4],v4;
                glGetMaterialfv(GL_FRONT, GL_AMBIENT,    v1);               // get default material
                glGetMaterialfv(GL_FRONT, GL_DIFFUSE,    v2);
                glGetMaterialfv(GL_FRONT, GL_SPECULAR,   v3);
                glGetMaterialfv(GL_FRONT, GL_SHININESS, &v4);

                float bndmat[4]= { bondmat[0]*4, bondmat[1]*4, bondmat[2]*4, 1.0 };
                glMaterialfv(GL_FRONT, GL_AMBIENT,  bndmat);                // set extra bond material
                glMaterialfv(GL_FRONT, GL_DIFFUSE,  bndmat);
                glMaterialfv(GL_FRONT, GL_SPECULAR, bndmat);
                glMaterialf (GL_FRONT, GL_SHININESS, 51.2f);

                glBegin(GL_LINES);
                for ( int i=0; i<bonds.size(); i++) {
                    if (select_mode == GL_SELECT) glLoadName(bn_level+i);
                    glNormal3d (ax1,ay1,az1);  glVertex3d (ax1,ay1,az1);
                    glNormal3d (ax2,ay2,az2);  glVertex3d (ax2,ay2,az2);
                }
                glEnd();

                glMaterialfv(GL_FRONT, GL_AMBIENT  , v1);                   // return back 
                glMaterialfv(GL_FRONT, GL_DIFFUSE  , v2);                   // default material
                glMaterialfv(GL_FRONT, GL_SPECULAR , v3);
                glMaterialf (GL_FRONT, GL_SHININESS, v4); 

            } else {                                                        // normal cylindrical bonds

                int i = 0;
                int sl = present ? slpresentcylinder : slnormalcylinder;

                for ( vector <struct bond>::iterator it  = bonds.begin(); 
                                                     it != bonds.end();    it++, i++) {
                    double o = (*it).o;
                    double l = len(ax1,ay1,az1,ax2,ay2,az2);
                    double a = 0.0;
                    if (l>0) a=acos((ay2-ay1)/l)/PI*180.0;
            
                    if(clrbnd) {                                        // twocolor bonds

                       glPushMatrix ();
                       if (select_mode == GL_SELECT) glLoadName(bn_level+i);
                       glTranslatef(ax1,ay1,az1);
                       glRotatef(a, az2-az1, 0.0, ax1-ax2);
                       glMaterialfv (GL_FRONT, GL_DIFFUSE, at_def[atoms[(*it).a-1].type-1]);
                       Cylinder(o,o,l/2,sl);
                       glPopMatrix ();

                       glPushMatrix ();
                       if (select_mode == GL_SELECT) glLoadName(bn_level+i);
                       glTranslatef((ax1+ax2)/2,(ay1+ay2)/2,(az1+az2)/2);
                       glRotatef(a, az2-az1, 0.0, ax1-ax2);
                       glMaterialfv (GL_FRONT, GL_DIFFUSE, at_def[atoms[(*it).b-1].type-1]);
                       Cylinder(o,o,l/2,sl);
                       glPopMatrix ();

                    } else {                                                // normal monocolor bonds

                       glPushMatrix ();
                       if (select_mode == GL_SELECT) glLoadName(bn_level+i);
                       glTranslatef(ax1,ay1,az1);
                       glRotatef(a, az2-az1, 0.0, ax1-ax2);
                       glMaterialfv(GL_FRONT, GL_DIFFUSE, bondmat);
                       Cylinder(o,o,l,sl);
                       glPopMatrix ();

                    }

                }   
          }
        }
        
        // Draw XYZ Axes
        if ( DrawAxes == 1 ) {
  
            float v2[4];
            glGetMaterialfv(GL_FRONT, GL_DIFFUSE,  v2);                    // get default material                   

            double l=maxcor*0.75;
            double d=maxcor*0.01;
            int sl = present ? slpresentcylinder : slnormalcylinder;

            glMaterialfv ( GL_FRONT, GL_DIFFUSE  , Black );                             // Origin
            glTranslatef(0.0, 0.0, 0.0);                                
            Sphere(d, sl, sl);

            glMaterialfv ( GL_FRONT, GL_DIFFUSE  , Red );                               // X axis
            glPushMatrix ();
            glTranslatef(0.0, 0.0, 0.0);                                                // origin
            glRotatef(270.0, 0, 0, 1);
            Cylinder( d, d, l, sl );
            glPopMatrix ();

            glPushMatrix ();
            glTranslatef(l, 0, 0);
            glRotatef(270.0, 0, 0, 1);
            Cylinder( d*2, 0.001,  d*6 , sl );
            Cylinder( d*2, d    , -d/5 , sl );
            glPopMatrix ();

    
            glMaterialfv ( GL_FRONT, GL_DIFFUSE  , Green );                             // Y axis
            glPushMatrix ();
            glTranslatef(0.0, 0.0, 0.0);                                                // origin
            Cylinder( d, d, l, sl );
            glPopMatrix ();

            glPushMatrix ();
            glTranslatef(0, l, 0);
            Cylinder( d*2, 0.001,  d*6 , sl );
            Cylinder( d*2, d    , -d/5 , sl );
            glPopMatrix ();
   

            glMaterialfv ( GL_FRONT, GL_DIFFUSE  , Blue );                              // Z axis
            glPushMatrix ();
            glTranslatef(0.0, 0.0, 0.0);                                                // origin
            glRotatef(90.0, 1, 0, 0);
            Cylinder( d, d, l, sl );
            glPopMatrix ();

            glPushMatrix ();
            glTranslatef(0, 0, l);
            glRotatef(90.0, 1, 0, 0);
            Cylinder( d*2, 0.001,  d*6 , sl );
            Cylinder( d*2, d    , -d/5 , sl );
            glPopMatrix ();

            glMaterialfv(GL_FRONT, GL_DIFFUSE  , v2);
        }

    	// Draw Gradient Vectors
        if(viewgr==1 && gradss==1){

            glMaterialfv(GL_FRONT, GL_DIFFUSE, irmat);

            int i;

            double vg=0;
            for(i=0;i<atoms.size();i++){                                   		// find largest gradient
              double l = len(0,0,0,grads[i].x,grads[i].y,grads[i].z);
              if(l>vg) vg=l;
            }
            vg=3/vg;

            for(i=0;i<atoms.size();i++) {

                  if (select_mode == GL_SELECT) glLoadName(gr_level+i);

                  double l = len(0,0,0,grads[i].x,grads[i].y,grads[i].z);

                  if (l>0) {                                                		// length of the vector shouldn't be zero

                       glPushMatrix ();
                       glTranslatef(atoms[i].x, atoms[i].y, atoms[i].z);
                       l*=vg;

                       double a=acos((vg*grads[i].y)/l)/PI*180.0;
                       glRotatef(a, vg*grads[i].z, 0, -vg*grads[i].x);

                       int sl = present ? slpresentcylinder : slnormalcylinder;

                       glPushMatrix();
                       glScaled(1,l,1);
                       Cylinder(0.05,0.05,1,sl);
                       glPopMatrix ();
                       glPopMatrix ();

                       glPushMatrix ();
                       glTranslatef(atoms[i].x+vg*grads[i].x, 
                                    atoms[i].y+vg*grads[i].y, 
                                    atoms[i].z+vg*grads[i].z);

                       glRotatef(a, vg*grads[i].z, 0, -vg*grads[i].x);

                       Cylinder( 0.1, 0.001,  0.3 , sl );				        // arrow
                       Cylinder( 0.1, 0.05 , -0.01, sl );
                       glPopMatrix ();

                  }
              }
        }

        glPopMatrix();

        glPopMatrix();
    }

    
    // HANDLE WINDOW RESIZING
    //    If window reshaped, need to readjust viewport/ortho
    void resize(int X,int Y,int ww,int hh) {
        if (ww>hh) ww=hh; else if (ww<hh) hh=ww;
        Fl_Gl_Window::resize(X,Y,ww,hh);
        glViewport (0, 0, ww, hh);
        SetOrtho();
        redraw();
    }

    
    int mouseX, mouseY, mouseS, mouseXo, mouseYo, mouseSo;
    #define MAXSELECT 256
    unsigned int *selectBuf;


    int handle (int event) {
        switch(event) {
            case FL_PUSH:
//              ... mouse down event ...
//              ... position in Fl::event_x() and Fl::event_y()
              if ( Fl::event_button() == 1) {

                  mouseX  = Fl::event_x(); 
                  mouseY  = Fl::event_y(); 
                  mouseS  = event;
                  cursor(FL_CURSOR_CROSS);

              }
              return 1;
            case FL_DRAG:
//              ... mouse moved while down event ...
              if ( Fl::event_button() == 1) {
                  mouseXo = mouseX;        mouseYo = mouseY;        mouseSo = mouseS;     
                  mouseX  = Fl::event_x(); mouseY  = Fl::event_y(); mouseS  = event;
                  rotxyz( 0,1,0, (mouseXo-mouseX),  roma );
                  rotxyz( 1,0,0, (mouseYo-mouseY),  roma );
                  redraw();
              }
              return 1;

            case FL_RELEASE: {   
//              ... mouse up event ...
                  cursor(FL_CURSOR_DEFAULT);
                                      
                  selectBuf = new unsigned int[MAXSELECT];
  
                  glSelectBuffer(MAXSELECT, selectBuf);
                  glRenderMode(GL_SELECT);
                  glInitNames();
                  glPushName((GLuint)~0);

                  x_select = mouseX;
                  y_select = mouseY;

                  select_mode = GL_SELECT;

                  glPushMatrix();

                  SetOrtho();
                  draw();

                  glPopMatrix();

                  int hits = glRenderMode(GL_RENDER); 

                  if (hits <= 0) {
                    select_hit = -1;
                  } else {
                    int o,i;
                    for( o=i=0; i<hits; i++ )
                      if (selectBuf[o*4+1] > selectBuf[i*4+1]) o=i;
                    select_hit = selectBuf[o*4+3];
                  }

                  delete[]selectBuf;

                  if (select_hit >= gr_level) {
                      int i = select_hit - gr_level;
                      double l = len(0,0,0,grads[i].x,grads[i].y,grads[i].z);
                      char s[64];
                      sprintf(s,"Gradient #%3d length is %12.6lf\n",i+1, l );
                      disp->insert_position( buff->length() );
                      disp->insert(s);
                  } else if (select_hit >= bn_level) {
                      int i = select_hit - bn_level;
                      ap1=bonds[i].a-1; ap2=bonds[i].b-1; 
                      char s[64];
                      sprintf(s,"Bond  (%3d,%3d) = %12.6lf\n",ap1+1,ap2+1, getBondLen(atoms[ap1],atoms[ap2]) );  
                      disp->insert_position( buff->length() );
                      disp->insert(s);
                  } else if (select_hit >= at_level) {
                      int i = select_hit - at_level;
                      ap4=ap3; ap3=ap2; ap2=ap1; ap1=i;
                      redraw();
                  }

                  select_mode = GL_RENDER;

                  SetOrtho();

              }
              return 1;

            case FL_MOUSEWHEEL :
              if (Fl::event_dy()<0) maxcor/=1.05;
              if (Fl::event_dy()>0) maxcor*=1.05;
              SetOrtho();
              redraw();
              zoom->value(maxcor);
              return 1;
            case FL_FOCUS :
            case FL_UNFOCUS :
//              ... Return 1 if you want keyboard events, 0 otherwise
              return 1;
            case FL_KEYBOARD:
//              ... keypress, key is in Fl::event_key(), ascii in Fl::event_text()
//              ... Return 1 if you understand/use the keyboard event, 0 otherwise...
              switch (Fl::event_key()) {
                 case FL_Escape :
                     win->hide();       // Quit
                 return 1;               
              }
              return 0;
            case FL_SHORTCUT:
//              ... shortcut, key is in Fl::event_key(), ascii in Fl::event_text()
//              ... Return 1 if you understand/use the shortcut event, 0 otherwise...
              return 0;
            default:
              // pass other events to the base class...
              return Fl_Gl_Window::handle(event);
            }
    }
    
    // Callback: when use picks 'File | Open' from main menu
    static void open_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->open_cb_real(w); }
    void open_cb_real(Fl_Widget*) {

        // Create the file chooser, and show it
        Fl_File_Chooser chooser(".",                                    // directory
        "XYZ Files (*.{xyz,mol})\tOUT Files (*.{out,log,gms,mop})",     // filter
        Fl_File_Chooser::SINGLE,                                        // chooser type
        "Open File");                                                   // title
        chooser.show();

        while(chooser.shown()) { Fl::wait(); }

        // User hit cancel?
        if ( chooser.value() != NULL ) { 

            atoms.clear();
            bonds.clear();
            freqs.clear();
            grads.clear();
            cread = 0;
            bread = 0;

            strcpy(fname, chooser.value());
            if ( ReadInputFile() ) {

                char title[256];
                strncpy(title, "ViewMol3D :: ", 256);
                strcat (title, fname);
                win->label(title);

                if ( bread == 0 ) CreateBond(mcalc,atoms_ehm);             // manual bonding
                CalcMaxCor();
                zoom->value(maxcor);
                InitGL();
                redraw();
            }

        }

    }


    // Callback: when user picks 'Quit'
    static void quit_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->quit_cb_real(w); }
    void quit_cb_real(Fl_Widget*) {
        win->hide();
    }


    // Callback: Perspective
    static void ren_pers_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->ren_pers_cb_real(w); }
    void ren_pers_cb_real(Fl_Widget*) {
        persp = 1 - persp;
        SetOrtho();
        redraw();
    }


    // Callback: Quality
    void RedrawQuality() {
        slnormalsphere    = int(32.*quality/10.);
        slnormalcylinder  = int(16.*quality/10.);
        slpresentsphere   = int(64.*quality/10.);
        slpresentcylinder = int(32.*quality/10.);
        redraw();
    }

    static void slider_q_callback(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->slider_q_callback_real(w); }
    void slider_q_callback_real(Fl_Widget* o) {
        quality = int(((Fl_Slider*)o)->value());
        RedrawQuality();
    }

    Fl_Window * winq;

    static void ok_butt_callback(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->ok_butt_callback_real(w); }
    void ok_butt_callback_real(Fl_Widget*o) {
        winq->hide();
    }

    static void cn_butt_callback(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->cn_butt_callback_real(w); }
    void cn_butt_callback_real(Fl_Widget*o) {
        quality = old_quality;              // restore quality state
        RedrawQuality();
        winq->hide();
    }

    static void ren_qual_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->ren_qual_cb_real(w); }
    void ren_qual_cb_real(Fl_Widget*) {

        int ww = 300, hh = 100;
        int bw = 80 , ds = 10;

        winq = new Fl_Window (ww, hh);
        winq->set_modal();

        Fl_Hor_Slider slider_q ( 0,  25,  ww,  30,  "Quality:");
        slider_q.align(FL_ALIGN_TOP);
        slider_q.callback(slider_q_callback,this);
        slider_q.value(quality);
        slider_q.step(1.0);
        slider_q.bounds(1.0, 20.0);

        Fl_Button ok_butt(1*ww/4-bw/2, 65, bw, 30, "Ok"    );  ok_butt.callback(ok_butt_callback, this);
        Fl_Button cn_butt(3*ww/4-bw/2, 65, bw, 30, "Cancel");  cn_butt.callback(cn_butt_callback, this);

        old_quality = quality;              // save quality state

        winq->show();

        Fl::run();

        winq->remove(slider_q);
        winq->remove(ok_butt);
        winq->remove(cn_butt);
        delete winq;
    }

    // Callback: Axes
    static void ren_axes_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->ren_axes_cb_real(w); }
    void ren_axes_cb_real(Fl_Widget*) { 
        DrawAxes = 1-DrawAxes;
    }

    // Callback: Reset Zoom
    static void ren_rszm_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->ren_rszm_cb_real(w); }
    void ren_rszm_cb_real(Fl_Widget*) {
        CalcMaxCor();
        SetOrtho();
        redraw();
    }

    // Callback: Gradients
    static void ren_grad_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->ren_grad_cb_real(w); }
    void ren_grad_cb_real(Fl_Widget*) {
        viewgr = 1 - viewgr;
        redraw();
    }

    
    // Callback: Standard Orientation
    static void ort_std_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->ort_std_cb_real(w); }
    void ort_std_cb_real(Fl_Widget*) {
        if ( atoms.size()>0 ) CalcPricipalAxes();
    }


    // Callback: Z-matrix Orientation
    static void ort_zmt_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->ort_zmt_cb_real(w); }
    void ort_zmt_cb_real(Fl_Widget*) {
        if ( atoms.size()>0 ) ZOrientMolecule( ap1, ap2, ap3 );
    }


    static void atm_off_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->atm_off_cb_real(w); }
    void atm_off_cb_real (Fl_Widget*) { bond  = 1-bond ; }

    static void atm_cpk_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->atm_cpk_cb_real(w); }
    void atm_cpk_cb_real (Fl_Widget*) { trad  = 1-trad ; }

    static void atm_num_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->atm_num_cb_real(w); }
    void atm_num_cb_real (Fl_Widget*) { atnum = 1-atnum; }

    static void atm_mar_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->atm_mar_cb_real(w); }
    void atm_mar_cb_real (Fl_Widget*) { mrkatm=1-mrkatm; }

    static void bnd_off_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->bnd_off_cb_real(w); }
    void bnd_off_cb_real (Fl_Widget*) { dline = 1-dline; }

    static void bnd_clr_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->bnd_clr_cb_real(w); }
    void bnd_clr_cb_real (Fl_Widget*) { clrbnd=1-clrbnd; }
    

    static void slider_x_callback(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->slider_x_callback_real(w); }
    void slider_x_callback_real(Fl_Widget* o) {
      angle_x = (((Fl_Slider*)o)->value());
      rotxyz(0,1,0, old_angle_x-angle_x, roma);
      old_angle_x=angle_x;
      redraw();
    }

    static void slider_y_callback(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->slider_y_callback_real(w); }
    void slider_y_callback_real(Fl_Widget* o) {
      angle_y = (((Fl_Slider*)o)->value());
      rotxyz(1,0,0, old_angle_y-angle_y, roma);
      old_angle_y=angle_y;
      redraw();
    }

    static void slider_z_callback(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->slider_z_callback_real(w); }
    void slider_z_callback_real(Fl_Widget* o) {
      angle_z = (((Fl_Slider*)o)->value());
      rotxyz(0,0,1, old_angle_z-angle_z, roma);
      old_angle_z=angle_z;
      redraw();
    }

    static void zoom_callback(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->zoom_callback_real(w); }
    void zoom_callback_real(Fl_Widget* o) {
      maxcor = (((Fl_Slider*)o)->value());
      SetOrtho();
      redraw();
    }

    static void about_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->about_cb_real(w); }
    void about_cb_real(Fl_Widget*) { 
        about->show();
    }

    static void getlen_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->getlen_cb_real(w); }
    void getlen_cb_real(Fl_Widget *) {
        char s[64];
        sprintf(s,"Length(%3u,%3u) = %12.6lf\n",ap1+1,ap2+1, 
                                            getBondLen(atoms[ap1],atoms[ap2]) );  
        disp->insert_position( buff->length() );
        disp->insert(s);
    }

    static void getang_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->getang_cb_real(w); }
    void getang_cb_real(Fl_Widget *) {
        char s[64];
        sprintf(s,"Angle (%3u,%3u,%3u) = %12.6lf\n",ap1+1,ap2+1,ap3+1, 
                                            getAngle(atoms[ap1],atoms[ap2],atoms[ap3]) );
        disp->insert_position( buff->length() );
        disp->insert(s);
    }

    static void getdih_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->getdih_cb_real(w); }
    void getdih_cb_real(Fl_Widget *){
        char s[64];
        sprintf(s,"Dihedr(%3u,%3u,%3u,%3u) = %12.6lf\n",ap1+1,ap2+1,ap3+1,ap4+1, 
                                            getDA(atoms[ap1],atoms[ap2],atoms[ap3],atoms[ap4]) );  
        disp->insert_position( buff->length() );
        disp->insert(s);
    }

    static void getzmt_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->getzmt_cb_real(w); }
    void getzmt_cb_real(Fl_Widget *) {
        char ss[4]; 
        ss[0]=at_nam[atoms[ap1].type-1][0];
        ss[1]=at_nam[atoms[ap1].type-1][1];
        if(ss[1]==' ') ss[1]=0; else ss[2]=0;
        double 
          bnd = getBondLen (atoms[ap1],atoms[ap2]),
          ang = getAngle   (atoms[ap1],atoms[ap2],atoms[ap3]),
          dih = getDA      (atoms[ap1],atoms[ap2],atoms[ap3],atoms[ap4]);
        char s[64];
        sprintf(s,"%s%02d  %3d %6.4lf  %3d %7.3lf  %3d %8.3lf\n",
                         ss,ap1+1, ap2+1,bnd, ap3+1,ang,  ap4+1,dih );
        disp->insert_position( buff->length() );
        disp->insert(s);
    }

    static void save_xyz_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->save_xyz_cb_real(w); }
    void save_xyz_cb_real(Fl_Widget *) {

        string namexyz = fname;
        namexyz += ".xyz";                                  // add to original file name XYZ extention

        ofstream fxyz( namexyz.c_str(), ios::out );         // open XYZ file for output

        fxyz << atoms.size() << endl;                       // number of atoms
        fxyz << "*** " << namexyz << endl;                  // comment: file name

        for ( vector <struct atom>::iterator atom  = atoms.begin();
                                             atom != atoms.end  ();  atom++ ) 
        {
/*            fxyz << setw(3)  << at_nam[(*atom).type-1]      // write atom name and its coorinates
                 << setw(16) << setprecision(7) << showpoint << (*atom).x 
                 << setw(16) << setprecision(7) << showpoint << (*atom).y 
                 << setw(16) << setprecision(7) << showpoint << (*atom).z 
                 << endl;*/                                 // C++ output style - looks pretty weird
                char s[256];
                sprintf(s,"%3s %16.7lf %16.7lf %16.7lf",
                        at_nam[(*atom).type-1], (*atom).x, (*atom).y, (*atom).z);
                fxyz << s << endl;
        }

        fxyz.close();                                       // that's it, close file

        namexyz = "XYZ file <" + namexyz + "> was saved.\n";
        disp->insert_position( buff->length() );
        disp->insert(namexyz.c_str());                      // print message

    }

    static void save_mol_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->save_mol_cb_real(w); }
    void save_mol_cb_real(Fl_Widget *) {

      char *namemol = new char[256];
      sprintf(namemol, "%s.mol", fname);

      FILE* fmol=fopen(namemol,"w");

      fprintf( fmol, "%3u ATOMS, %3u BONDS,    0 CHARGES, %s\n", atoms.size(), bonds.size(), fname);

      vector <struct atom>::iterator atom;
      int i = 1;
      for ( atom = atoms.begin(); atom!=atoms.end(); atom++, i++ ) {
         fprintf( fmol, "%5u %3s %12.8lf %12.8lf %12.8lf    0.000\n", i, at_nam[(*atom).type-1], (*atom).x, (*atom).y, (*atom).z);
      }

      vector <struct bond>::iterator bond;
          i = 1;
      for ( bond = bonds.begin(); bond!=bonds.end(); bond++, i++ ) {
         fprintf( fmol, "%5u  %6u  %6u     SINGLE\n", i, (*bond).a, (*bond).b);
      }

      fclose(fmol);

      char msg[256];
      sprintf(msg,"MOL file <%s> was saved.\n", namemol);
      disp->insert_position( buff->length() );
      disp->insert(msg);                  // print message

    }

    // MINDO/3 block, energy, gradients calculation and optimization
    static void mindo3_eng_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->mindo3_eng_cb_real(w); }
    void mindo3_eng_cb_real(Fl_Widget *w) {

        struct atom_ehm at;

        atoms_ehm.erase( atoms_ehm.begin(), atoms_ehm.end() );
        for (vector <struct atom>::iterator atit = atoms.begin(); atit!=atoms.end(); atit++) {
          at.atno =   (*atit).type;
          at.xyz[0] = (*atit).x;
          at.xyz[1] = (*atit).y;
          at.xyz[2] = (*atit).z;
          atoms_ehm.push_back(at);
        }

        if (mcalc != NULL) delete mcalc;

        mcalc = new mindo3(&atoms_ehm);

    //  ehm_Eng = mcalc->ehm::Energy();
        Energy = mcalc->Energy();

        char msg[256];
        sprintf(msg,"SCF (%u/%u) done. Ediff=%12.6le\nHeat of Formation = %9.4f kcal/mol\n", mcalc->SCFit, mcalc->MaxSCF, mcalc->Ediff, Energy);
        disp->insert_position( buff->length() );
        disp->insert(msg);                  // print message

    }

    static void mindo3_grd_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->mindo3_grd_cb_real(w); }
    void mindo3_grd_cb_real(Fl_Widget *w) {

        Matrix G;

        mindo3_eng_cb_real(w);
        mcalc->forces(&G);            //  mcalc->num_forces(atoms_ehm, &G,0.00000001);
      //  mcalc->num_forces_central4(&G,1e-8);

        struct grad grad;
        int i;

        for ( i=0; i<atoms.size(); i++) {
          grad.znuc =  atoms[i].type;
          grad.x    = -G(i+1,1); 
          grad.y    = -G(i+1,2); 
          grad.z    = -G(i+1,3); 
          if (gradss == 1) grads[i] = grad;
          else             grads.push_back( grad );
        }  

        gradss = 1;

        int nmaxg=1;
        double maxg=0.0, ssum=0.0;
        for( i=0; i<atoms.size(); i++ ){
          double l = len(0,0,0, grads[i].x,grads[i].y,grads[i].z);
          if (l>maxg) { maxg=l; nmaxg=i; };
          ssum += sqr(G(i+1,1))+sqr(G(i+1,2))+sqr(G(i+1,3));
        }
        RMS = sqrt(ssum/(atoms.size()*3));

        char msg[256];
        sprintf(msg,"RMS = %9.6lf, Maximum gradient on atom %3d is\ndE/dx=%8.5f dE/dy=%8.5f dE/dz=%8.5f\n",
                     RMS, nmaxg+1, grads[nmaxg].x, grads[nmaxg].y, grads[nmaxg].z);
        disp->insert_position( buff->length() );
        disp->insert(msg);                  // print message

    }

    Fl_Window           *win_progress;
    Fl_Progress         *progress;
    Fl_Multiline_Output *output;
    Fl_Button           *prcn_butt; 

    int CancelOptimization;

    static void prcn_butt_callback(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->prcn_butt_callback_real(w); }
    void prcn_butt_callback_real(Fl_Widget * w) { 
        CancelOptimization = TRUE;
    }
            
    void funcgrad(ap::real_1d_array x, double& f, ap::real_1d_array& g) {

        Matrix G;
        int i;

        if (CancelOptimization) {
            f = 0.0;
            for ( i=0; i<atoms.size(); i++ ) {
              g(i*3+1) = 0,
              g(i*3+2) = 0,
              g(i*3+3) = 0;
            }
            return;
        }

      for (i=0; i<atoms.size(); i++) {
        atoms_ehm[i].xyz[0] = x(i*3+1);
        atoms_ehm[i].xyz[1] = x(i*3+2);
        atoms_ehm[i].xyz[2] = x(i*3+3);
      }

      mcalc->D_initialized  = FALSE;
      mcalc->F0_initialized = FALSE;

      Energy = mcalc->Energy(); 
      f = Energy;
      mcalc->forces(&G);

      double ssum = 0.0;
      for ( i=0; i<atoms.size(); i++ ) {
        g(i*3+1) = G(i+1,1),
        g(i*3+2) = G(i+1,2),
        g(i*3+3) = G(i+1,3);
        ssum += sqr(G(i+1,1))+sqr(G(i+1,2))+sqr(G(i+1,3));
      }
      RMS = sqrt(ssum/(atoms.size()*3));
  
      funcalls++;

      char msg[256];
      sprintf(msg,"Energy = %9.6f kcal/mol\nafter SCF (%2u/%2u)\nRMS[%3u] = %9.6lf", 
                   f,               mcalc->SCFit, mcalc->MaxSCF, funcalls,  RMS);
      output->value(msg);

      progress->value(double(funcalls)/32.0);
      win_progress->redraw();

      for ( i=0; i<atoms.size(); i++ ) {
        atoms[i].x = x(i*3+1);
        atoms[i].y = x(i*3+2);
        atoms[i].z = x(i*3+3);
      }

      for ( i=0; i<atoms.size(); i++) {
        struct grad grad;
        grad.znuc =  atoms[i].type;
        grad.x    = -G(i+1,1); 
        grad.y    = -G(i+1,2); 
        grad.z    = -G(i+1,3); 
        if (gradss == 1) grads[i] = grad;
        else             grads.push_back( grad );
      }  
      gradss = 1;

      redraw();

      Fl::check();

    }

    static void mindo3_opt_cb(Fl_Widget* w, void* data) { ((MyGlWindow*)data)->mindo3_opt_cb_real(w); }
    void mindo3_opt_cb_real(Fl_Widget * w) {
      struct atom_ehm at;
      int i, n = atoms.size()*3, m = 5;
      ap::real_1d_array x;
      double eps = 0.0001;

      const int pw = 250, ph = 80, pp = 30, oo = 55, bb = 25;

      win_progress = new Fl_Window(pw,120);
      win_progress->set_modal();
      progress = new Fl_Progress(0,0,pw,pp);
      progress->minimum(0);               // set progress bar attribs..
      progress->maximum(1);
      output = new Fl_Multiline_Output(0,pp+5,pw,oo); 
      output->value("Start optimization...");
      CancelOptimization = FALSE;
      prcn_butt = new Fl_Button(0, pp+oo+10, pw, bb, "Cancel");  
      prcn_butt->callback(prcn_butt_callback,this);

      win_progress->show();
      win_progress->redraw();
      Fl::check();

      x.setbounds(1, n);

      CenterMolecule();

      atoms_ehm.erase( atoms_ehm.begin(), atoms_ehm.end() );
      for ( i=0; i<atoms.size(); i++ ) {
        x(i*3+1)  = atoms[i].x,
        x(i*3+2)  = atoms[i].y,
        x(i*3+3)  = atoms[i].z;
        at.atno   = atoms[i].type;
        at.xyz[0] = atoms[i].x;
        at.xyz[1] = atoms[i].y;
        at.xyz[2] = atoms[i].z;
        atoms_ehm.push_back(at);
      }

      if (mcalc != NULL) delete mcalc;
      mcalc = new mindo3(&atoms_ehm);

      funcalls = 0;

      int result = lbfgsminimize(n, m, x, eps);

      char msg[256];  
      if (CancelOptimization) {
        sprintf(msg,"Optimization was canceled, %u iterations\n",funcalls);
        disp->insert_position( buff->length() );
        disp->insert(msg);                  // print message
      } else if (result==0) {
        sprintf(msg,"Optimization completed successfully, %u iterations\n",funcalls);
        disp->insert_position( buff->length() );
        disp->insert(msg);                  // print message
      } else {
        sprintf(msg,"Optimization failed with error: %u, %u iterations\n",result,funcalls);
        disp->insert_position( buff->length() );
        disp->insert(msg);                  // print message
      }
      sprintf(msg,"Final energy = %12.6f kcal/mol, RMS = %9.6lf\n", Energy, RMS);
      disp->insert_position( buff->length() );
      disp->insert(msg);                  // print message

      for ( i=0; i<atoms.size(); i++ ) {
        atoms[i].x = x(i*3+1);
        atoms[i].y = x(i*3+2);
        atoms[i].z = x(i*3+3);
      }

      win_progress->remove(output);                 // remove output
      win_progress->remove(progress);               // remove progress bar from window
      win_progress->remove(prcn_butt);                // remove 
      delete output;
      delete progress;                   		    // deallocate it
      delete prcn_butt;
      delete win_progress;

    }

    // Constructor
    MyGlWindow(Fl_Window* _win, char *_fname, int X, int Y, int W, int H, const char* L=0) : Fl_Gl_Window(X,Y,W,H,L),
        trad(0), bond(0), persp(0), LowQual(0), dline(0), chg(0), present(0), monum(1), atnum(0), viewir(0),      
        frcur(0), viewgr(0), strnum(0), xaccl(0), screenNum(1), irc(0), clrbnd(0), pbuff(0), mrkatm(0),
        DrawMO(0), DrawAxes(0), DrawSurf(0), ap1(0), ap2(1), ap3(2), ap4(3), quality(10),
        slnormalsphere(32), slnormalcylinder(16), slpresentsphere(64), slpresentcylinder(32),
        select_mode(GL_RENDER), select_hit(0), numorb(0),
        old_angle_x(0), angle_x(0), old_angle_y(0), angle_y(0), old_angle_z(0), angle_z(0), win(_win)
    {

        win->resizable(this);
        mode(FL_RGB | FL_DOUBLE | FL_DEPTH);

        strncpy(fname, _fname, 256);

        int i,j;

        for (i=strlen(fname)-1;i!=0;i--)                    // delete trailing spaces
            if (fname[i]==' ') fname[i]=0; 
            else break;

#ifdef WIN32
        char tmp[256];                                          // delete both '"'
        j=0;
        for ( int k=0; k<strlen(fname); k++ )
            if ( fname[k]!='"' ) { tmp[j]=fname[k]; j++; }
        tmp[j] = 0;
        strncpy( fname, tmp, 256);
#endif

        if ( strlen(fname) != 0 ) {
            if ( ReadInputFile() ) {
                if ( bread == 0 ) CreateBond(mcalc,atoms_ehm);
                CalcMaxCor();
            }
        }

        for (i=0; i<4; i++)
            for (j=0; j<4; j++) 
                if ( i==j) roma[i][j] = 1; else roma[i][j] = 0;

        menubar = new Fl_Menu_Bar (0,0, win->w()-110,sw);

        menubar->add("&File/&Open",     0, open_cb    , this);
        menubar->add("&File/Save &XYZ", 0, save_xyz_cb, this);
        menubar->add("&File/Save &MOL", 0, save_mol_cb, this);
        menubar->add("&File/&Quit",     0, quit_cb    , this);

        menubar->add("&Render/&Perspective", 0, ren_pers_cb, this, FL_MENU_TOGGLE);
        menubar->add("&Render/&Quality"    , 0, ren_qual_cb, this);
        menubar->add("&Render/&Axes"       , 0, ren_axes_cb, this, FL_MENU_TOGGLE);
        menubar->add("&Render/&Reset zoom" , 0, ren_rszm_cb, this, 0);
        menubar->add("&Render/&Gradients"  , 0, ren_grad_cb, this, FL_MENU_TOGGLE);

        menubar->add("&Atoms/&Off",     0, atm_off_cb, this, FL_MENU_TOGGLE);
        menubar->add("&Atoms/&CPK",     0, atm_cpk_cb, this, FL_MENU_TOGGLE);
        menubar->add("&Atoms/&Numbers", 0, atm_num_cb, this, FL_MENU_TOGGLE);
        menubar->add("&Atoms/&Marked",  0, atm_mar_cb, this, FL_MENU_TOGGLE);

        menubar->add("&Bonds/&Off"  , 0, bnd_off_cb, this, FL_MENU_TOGGLE);
        menubar->add("&Bonds/&Color", 0, bnd_clr_cb, this, FL_MENU_TOGGLE);

        menubar->add("&Orientation/&Standard", 0, ort_std_cb, this);
        menubar->add("&Orientation/&Z-Matrix", 0, ort_zmt_cb, this);

        menubar->add("&MINDO3/&Energy"      , 0, mindo3_eng_cb, this);
        menubar->add("&MINDO3/&Gradients"   , 0, mindo3_grd_cb, this);
        menubar->add("&MINDO3/&Optimization", 0, mindo3_opt_cb, this);

        menubar->add("A&bout", 0, about_cb, this);

        slider_x = new Fl_Hor_Slider  ( sw-ds,  win->h()-sw-2*ds,  win->w()-sw-ds-tw,   sw+ds,  "X:");
        slider_x->align(FL_ALIGN_LEFT);
        slider_x->callback(slider_x_callback, this);
        slider_x->value(angle_x);
        slider_x->step(1);
        slider_x->bounds(0,360);

        slider_y = new Fl_Slider (  0,  sw+3*ds,  sw+ds,  win->h()-(sw+3*ds)-(sw+3*ds),  "Y:");
        slider_y->align(FL_ALIGN_TOP);
        slider_y->callback(slider_y_callback, this);
        slider_y->value(angle_y);
        slider_y->step(1);
        slider_y->bounds(0,360);

        slider_z = new Fl_Slider ( win->w()-(sw+ds)-tw,   sw+3*ds,   sw+ds,   win->h()-(sw+3*ds)-(sw+3*ds),  "Z:");
        slider_z->align(FL_ALIGN_TOP);
        slider_z->callback(slider_z_callback, this);
        slider_z->value(angle_z);
        slider_z->step(1);
        slider_z->bounds(0,360);

        zoom = new Fl_Roller ( win->w()-rw,  0,  rw,  sr+ds,  "Zoom:");
        zoom->align(FL_ALIGN_LEFT);
        zoom->callback(zoom_callback, this);
        zoom->value(maxcor);
        zoom->step(0.1);
        zoom->minimum(0);
        zoom->maximum(100);

        // Information message window
        buff = new Fl_Text_Buffer();
        disp = new Fl_Text_Display(win->w()-tw+ds, 3*sw+4*ds, tw-ds, win->h()-(3*sw+5*ds), "Information:");
        disp->textfont(FL_COURIER);
        disp->textsize(12);
        disp->buffer(buff);
        //win->resizable(disp);

        // Buttons
        int bx0 = win->w()-tw+ds,                   // x0 for first button
            bw  = tw/4-ds;                          // button width
        get_len = new Fl_Button(bx0          , sw+4*ds, bw, sw, "&Length");
        get_ang = new Fl_Button(bx0+  bw+  ds, sw+4*ds, bw, sw, "&Angle");
        get_dih = new Fl_Button(bx0+2*bw+2*ds, sw+4*ds, bw, sw, "&Dihedral");
        get_zmt = new Fl_Button(bx0+3*bw+3*ds, sw+4*ds, bw, sw, "&Z-matrix");
        get_len->callback(getlen_cb, this);
        get_ang->callback(getang_cb, this);
        get_dih->callback(getdih_cb, this);
        get_zmt->callback(getzmt_cb, this);

        about = new Fl_Help_Dialog();
        about->value(AboutText);
    }

    // Destructor
    ~MyGlWindow() 
    {
        win->remove(menubar);           delete menubar;
        win->remove(slider_x);          delete slider_x;
        win->remove(slider_y);          delete slider_y;
        win->remove(slider_z);          delete slider_z;
        win->remove(zoom);              delete zoom;
        win->remove(get_len);           delete get_len;
        win->remove(get_ang);           delete get_ang;
        win->remove(get_dih);           delete get_dih;
        win->remove(get_zmt);           delete get_zmt;
        delete disp;
        delete buff;
        delete about;
    }

};


#ifdef WIN32
int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpszCmdParam, int nCmdShow) {
    char fname[256];
    strncpy(fname, lpszCmdParam, 256);
#else
int main(int argc, char* argv[]) {
    char fname[256];
    if (argc>1) strncpy(fname, argv[1], 256); else strncpy(fname, "", 256);
#endif

    char title[256];
    strncpy(title, "ViewMol3D :: ", 256);
    strcat (title, fname);

    // create main window
    Fl_Window  win  ( ww, wh, title );

    // create OpenGL window   
    MyGlWindow  vm3w ( &win, fname, (sw+3*ds), (sw+3*ds), win.w()-2*(sw+3*ds)-tw, win.h()-2*(sw+3*ds) );

    win.show();

    int result = Fl::run();
    
    return result;

}  
