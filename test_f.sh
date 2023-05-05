gcc -Wall -I/usr/local/include -O2 -c lua_fsolve.c
gcc lua_fsolve.o -O2 -shared -o params.dll -fPIC -llua -lgsl -lgslcblas -lm
rm lua_fsolve.o
