.PHONY: all

all:
	gfortran -Ofast -march=native -fPIC -shared -o potri2.so dpotri2.f90 -L/usr/local/lib -llapack -lblas -lgfortran -lm

mkl:
	gfortran -Ofast -march=native -fPIC -shared -o potri2.so dpotri2.f90 -L/usr/lib/x86_64-linux-gnu/mkl -llapack -lblas -lgfortran -lm
