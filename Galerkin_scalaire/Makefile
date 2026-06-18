FC	= gfortran
LD	= gfortran

FFLAGS  = -g -C -O3 -march=native 
LDFLAGS = -g -C -O3 -march=native 


# FFLAGS  = -g -C -O3 -fbacktrace -Wtarget-lifetime -Wcompare-reals -Wrealloc-lhs-all \
# 			-Wrealloc-lhs -Wintrinsic-shadow -Wunderflow -Wintrinsics-std \
# 			-Wimplicit-procedure -Wconversion -Wconversion-extra -Wline-truncation \
# 			-Wcharacter-truncation -Wc-binding-type -Waliasing -ffpe-trap=zero,overflow,underflow
# LDFLAGS = -g -C -O3 -fbacktrace -Wtarget-lifetime -Wcompare-reals -Wrealloc-lhs-all \
# 			-Wrealloc-lhs -Wintrinsic-shadow -Wunderflow -Wintrinsics-std \
# 			-Wimplicit-procedure -Wconversion -Wconversion-extra -Wline-truncation \
# 			-Wcharacter-truncation -Wc-binding-type -Waliasing -ffpe-trap=zero,overflow,underflow 

LIB = 
LIBPATH = 
MODPATH =

PROG = DG.exe

SRCS =	precis.F90 mod_declaration.F90 mod_init.F90 mod_Monolith.F90 \
	    mod_RKDG.F90 mod_polynome.F90 tests.F90 main.F90 mod_divers.F90

OBJS =	precis.o mod_declaration.o mod_init.o mod_Monolith.o \
		mod_RKDG.o mod_polynome.o main.o mod_divers.o 

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS)  $(LIBPATH)  $(MODPATH) -o $@ $(OBJS) $(LIB)
clean:
	rm -f $(PROG) $(OBJS) *.mod core *~ .*~ \#*# .#*

.SUFFIXES: .F90

%.o : %.F90
	$(FC) $(FFLAGS) $(MODPATH) -c $*.F90

mod_declaration.o : precis.o   
mod_polynome.o : mod_declaration.o
mod_divers.o : mod_declaration.o 
mod_Monolith.o : mod_polynome.o
mod_RKDG.o : mod_Monolith.o mod_divers.o 
mod_init.o : mod_polynome.o 
main.o : mod_RKDG.o mod_init.o 
# tests.o : mod_RKDG.o mod_init.o