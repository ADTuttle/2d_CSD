
PETSC_DIR=/Users/Austin/Documents/C_Progs/petsc-3.7.7/arch-darwin-c-debug

# include ${PETSC_DIR}/lib/petsc/conf/variables
# include ${PETSC_DIR}/lib/petsc/conf/rules

Includes = -I$(PETSC_DIR)/../include -I$(PETSC_DIR)/include -I/opt/X11/include -I$(PETSC_DIR)/../include/petsc/mpiuni -I$(PETSC_DIR)/include/petscconf.h

CSD_Includes = -I constants.h -I functions.h

Flags = -Wall -march=native -mtune=native

Linker= -L$(PETSC_DIR)/lib -L/opt/X11/lib

LFlags= -Wall -march=native -mtune=native 
csd: 
	gcc -O3 -o csd_main.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/csd_main.c
	gcc -O3 -o initialize.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/initialize.c
	gcc -O3 -o ion_channel.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/ion_channel.c
	gcc -O3 -o array_function.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/array_functions.c
	gcc -O3 -o update_solution.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/update_solution.c


	gcc -O3 $(LFlags)   -o csd initialize.o csd_main.o ion_channel.o array_function.o update_solution.o $(Linker) -lpetsc -lf2clapack -lf2cblas -lX11 -ldl 
	#View assembly file?
	# gcc -O3 $(LFlags)   -S csd initialize.o csd_main.o ion_channel.o array_function.o update_solution.o $(Linker) -lpetsc -lf2clapack -lf2cblas -lX11 -ldl 

	rm csd_main.o initialize.o ion_channel.o array_function.o update_solution.o

debug: 
	clang -g -o csd_main.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/csd_main.c
	clang -g -o initialize.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/initialize.c
	clang -g -o ion_channel.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/ion_channel.c
	clang -g -o array_function.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/array_functions.c
	clang -g -o update_solution.o -c $(Flags) $(Includes) $(CSD_Includes) `pwd`/update_solution.c


	clang -g $(LFlags)   -o csd initialize.o csd_main.o ion_channel.o array_function.o update_solution.o $(Linker) -lpetsc -lf2clapack -lf2cblas -lX11 -ldl 
	#View assembly file?
	# gcc -O3 $(LFlags)   -S csd initialize.o csd_main.o ion_channel.o array_function.o update_solution.o $(Linker) -lpetsc -lf2clapack -lf2cblas -lX11 -ldl 

	rm csd_main.o initialize.o ion_channel.o array_function.o update_solution.o

ex1: ex1.o  chkopts
	-${CLINKER} -o ex1 ex1.o  ${PETSC_KSP_LIB}
	${RM} ex1.o


