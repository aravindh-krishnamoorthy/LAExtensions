all: dtrgdsy dtrtrsy

dtrgdsy: dtrgdsy.so dtrgdsy.asm
dtrtrsy: dtrtrsy.so dtrtrsy.asm

clean:
	rm -rf dtrgdsy.so dtrtrsy.so

%.so: %.f90
	gfortran -Ofast -shared -march=native -o $@ $<

%.asm: %.f90
	gfortran -Ofast -S -masm=intel -march=native -o $@ $<
