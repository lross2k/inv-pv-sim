gcc -Wall -I/usr/local/include -c lua_fdsolve.c
gcc lua_fdsolve.o -shared -o newlib.dll -fPIC -llua -lgsl -lgslcblas -lm
