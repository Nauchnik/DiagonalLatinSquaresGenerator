CPP = mpiicpc
CPPFLAGS = -O3 -std=c++0x -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D _MPI

odls_finder: odls_finder.o odls_sequential.o
	${CPP} ${CPPFLAGS} odls_finder.o odls_sequential.o -o odls_finder
	
odls_finder.o: odls_finder.cpp
	${CPP} ${CPPFLAGS} odls_finder.cpp -c

odls_sequential.o: odls_sequential.cpp
	${CPP} ${CPPFLAGS} odls_sequential.cpp -c
	
clean:
	rm -rf *.o
	rm odls_finder
	clear