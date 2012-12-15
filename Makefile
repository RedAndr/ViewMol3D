# ViewMol3D
# Makefile for GCC compiler

CPP  = gcc

INCLUDE= -Inewmat 
CFLAGS = -c -O2
RFLAGS = -o vm3

.SUFFIXES : .o .cpp

.cpp.o :
	$(CPP) $(CFLAGS) $(INCLUDE) $<

OBJS   = cr_bond.o  ehm.o  geomath.o  graphprim.o  main.o  mindo3.o  readfile.o    
LFLAGS = -s -Lnewmat -L/usr/X11R6/lib/ -L/lib/ -static
MYLIBS = -lnewmat -lfltk_images -lfltk_gl -lfltk
LIBS   = -lpthread -lpng -ljpeg -lz -lglut -lGLU -lGL -ldl -lXext -lXxf86vm  -lX11 -lm -lc -lgcc -lstdc++ 

vm3: $(OBJS)
	$(CPP) $(RFLAGS) $(OBJS) $(LFLAGS) $(MYLIBS) $(LIBS) 
#$@ $?

main.o                  :  atom_rads.h debug.h gl_param.h units.h   main.cpp
cr_bond.o               :  atom_rads.h debug.h gl_param.h units.h   cr_bond.h cr_bond.cpp
ehm.o                   :  atom_rads.h debug.h gl_param.h units.h   ehm.h ehm.cpp
geomath.o               :  atom_rads.h debug.h gl_param.h units.h   geomath.h geomath.cpp
graphprim.o             :  atom_rads.h debug.h gl_param.h units.h   graphprim.h graphprim.cpp 
mindo3.o                :  atom_rads.h debug.h gl_param.h units.h   mindo3.h mindo3.cpp
readfile.o              :  atom_rads.h debug.h gl_param.h units.h   readfile.h readfile.cpp

all: vm3

clean:
	rm *.o *~ *.core vm3
