BASE = $(shell /bin/pwd)
EXED = $(BASE)/$(EXECTBL)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src
INCD = $(BASE)/header

#Compiler of C
CC=gcc

#Flags -I tells which directory to look for the .h files 
CFLAGS =-Wall -I$(INCD) -O3 $(shell gsl-config --cflags)

# GSL Libraries
LDFLAGS = $(shell gsl-config --libs)

#Libraries (Math)
LIBS=

#set of headers .h files
DEPS = $(INCD)/$(wildcard *.h)


#set of source .c files
SRC = load_data.c utilities.c debug.c spectral_utilities.c Butcher_Tableau.c Get_Differentiation_Matrices.c derivatives.c allocate_and_free_derivs.c solve_equations.c parameters.c hyperboloidal_functions.c Coordinates.c IndexRoutines.c Get_Initial_Guess.c newton_SDIRK_IG.c 1+1_Free_Boundaries.c FieldAndBoundEqns.c FuncAndJacobian.c newton.c newton_direct.c SDIRK_Preconditioner.c output.c get_Initial_Data.c
#set of object .o files
OBJS = $(SRC:.c=.o)

EXECTBL  = Wave

#rule to say that .o files depends on .c and files on DEPS. All those
#files are going to be checked on the execution for the up-to-date version


all:$(EXECTBL)


$(EXECTBL): $(OBJS)
	@echo '' 
	@echo '============ ...building the executable ============'
	@echo ''
	cd $(OBJD); $(CC) $(CFLAGS) -o ../$(EXECTBL) $(OBJS) $(LIBS) $(LDFLAGS)
	@echo ''
	@echo '============ ...Done! ============'


vpath %.o $(OBJD)
vpath %.c $(SRCD)
%.o:%.c $(DEPS)
	@echo '' 
	@echo '============ ...building' $*.o '============'
	@echo ''
	cd $(OBJD); $(CC) -o $@ -c $(SRCD)/$(*).c $(LIBS) $(CFLAGS)



flush:
	rm -f $(EXECTBL) *.dat fort.* $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.d $(OBJD)/*.pc*

clean:
	rm -f *.dat *.yg fort.* $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.d $(OBJD)/*.pc* $(SRCD)/*.o $(EXECTBL)
	rm -f *~ $(OBJD)/*~ $(SRCD)/*~ header/*~