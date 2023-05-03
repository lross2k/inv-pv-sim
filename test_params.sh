gcc -Wall -I/usr/local/include -c lua_fsolve.c
gcc lua_fsolve.o -shared -o newlib.dll -fPIC -llua -lgsl -lgslcblas -lm
