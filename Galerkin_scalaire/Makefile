#
FC=gfortran -g -C -O3 -march=native

DG : main.o  precis.o mod_declaration.o mod_polynome.o mod_RKDG.o mod_init.o
	$(FC) -o DG.exe main.o mod_RKDG.o mod_init.o mod_declaration.o mod_polynome.o precis.o  


tests : tests.o  precis.o mod_declaration.o mod_polynome.o mod_RKDG.o mod_init.o
	$(FC) -o tests.exe tests.o mod_RKDG.o mod_init.o mod_declaration.o mod_polynome.o precis.o  

precis.o : precis.F90
	$(FC) -c precis.F90

mod_declaration.o : precis.o mod_declaration.F90
	$(FC) -c mod_declaration.F90

mod_polynome.o : mod_declaration.o  mod_polynome.F90
	$(FC) -c mod_polynome.F90
	
mod_RKDG.o : mod_polynome.o  mod_RKDG.F90
	$(FC) -c mod_RKDG.F90
	
mod_init.o : mod_polynome.o  mod_init.F90
	$(FC) -c mod_init.F90

main.o : mod_RKDG.o mod_init.o main.F90 
	$(FC) -c main.F90
	
tests.o : mod_RKDG.o mod_init.o tests.F90 
	$(FC) -c tests.F90

clear : 
	rm -f *.mod 
	rm -f *.o
	rm -f *.exe