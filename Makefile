CC = gcc
DLLFLAGS = -shared -fPIC
CFLAGS = -Wall -Werror -O2 -I/mingw64/include
LINKS = -L/mingw64/lib -llua -lgsl -lgslcblas -lm

quick_test:
	$(CC) $(CFLAGS) -c call_from_lua.c
	$(CC) call_from_lua.o $(CFLAGS) $(DLLFLAGS) -L/mingw64/lib -llua -o NEWLIB.DLL

params_based:
	$(CC) lua_fdfsolve.c $(CFLAGS) -c -o fdf_solver.o
	#$(CC) fdf_solver.o $(CFLAGS) $(DLLFLAGS) $(LINKS) -o NEWLIB.DLL
	$(CC) fdf_solver.o $(CFLAGS) -bundle -undefined dynamic_lookup -fPIC $(LINKS) -o NEWLIB.DLL

two_diode:
	$(CC) $(CFLAGS) $(LINKS) -c -o two_diode_solver.o two_diode_solver.c
	$(CC) two_diode_solver.o $(CFLAGS) $(DLLFLAGS) $(LINKS) -o TDS.DLL 

