all: dtrgdsy.so dtrtrsy.so

clean:
	rm -rf dtrgdsy.so dtrtrsy.so

%.so: %.f90
	gfortran -Ofast -shared -masm=intel -march=native -o $@ $<
