gcc -Wall -I/usr/local/include -O2 -c lua_fdfsolve.c
gcc lua_fdfsolve.o -O2 -shared -o params.dll -fPIC -llua -lgsl -lgslcblas -lm
rm lua_fdfsolve.o
