.PHONY: all

SHLIB := potri2.so
SRCS  := dpotri2.f90 dpotri2s.f90
FOPTS += -Ofast -march=native -masm=intel -fverbose-asm -fPIC -shared -save-temps
FLIBS := -llapack -lblas -lgfortran -lm

all:
	gfortran $(FOPTS) -o $(SHLIB) $(SRCS) -L/usr/local/lib $(FLIBS)

mkl:
	gfortran $(FOPTS) -o $(SHLIB) $(SRCS) -L/usr/lib/x86_64-linux-gnu/mkl $(FLIBS)
