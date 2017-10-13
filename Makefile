SOURCES= lbm.cpp domain.cpp solutionExporterEulerian.cpp lattice.cpp function.cpp geometry.cpp boundaryc.cpp compute_quantities.cpp relaxation_time.cpp rarefied_models.cpp

OBJECTS=$(SOURCES:.cpp=.o)

CC=g++
LIBS=
LIB_PATHS=
INCLUDE_PATHS=
CFLAGS=-O2 -std=c++11 -Wall	 -W -Wno-uninitialized -fstrict-aliasing

EXECUTABLE=lbm

MAKE_CMD=$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIB_PATHS) $(INCLUDE_PATHS) $(LIBS) 

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(MAKE_CMD)

.cpp.o:
	$(CC) $(CFLAGS) -c -o $@ $<

remake:
	rm $(EXECUTABLE)
	$(MAKE_CMD)

clean:
	rm $(EXECUTABLE)
	rm $(OBJECTS)

omp: lbm.o  domain.o solutionExporterEulerian.o lattice.o function.o geometry.o boundaryc.o compute_quantities.o relaxation_time.o rarefied_models.o
	$(CC) $(CFLAGS) -o lbm lbm.o  domain.o solutionExporterEulerian.o lattice.o function.o geometry.o boundaryc.o compute_quantities.o relaxation_time.o rarefied_models.o $(LIB_PATHS) $(INCLUDE_PATHS) $(LIBS) 


