MODULE=pyrazinef90

gnu: FCOMPILER=gnu95
gnu: f2py

ifort: FCOMPILER=intelem
ifort: f2py2

f2py: $(MODULE).so

$(MODULE).so: pot_24.f90
	f2py --fcompiler=$(FCOMPILER) -c -m $(MODULE) $^

clean:
	$(RM) -f *.o *.so *.mod *.c *.pyc
