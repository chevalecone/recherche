CC = g++
CFLAGS = -std=c++11 -O2 -fstrict-aliasing -Wall -W -Wno-uninitialized 
EXEC_NAME = lbm
INCLUDES = 
LIBS = libtinyxml2.a
OBJ_FILES = lbm.o  domain.o parser.o tinyxml2.o solutionExporterEulerian.o lattice.o function.o geometry.o boundaryc.o compute_quantities.o relaxation_time.o

all : staticlib $(EXEC_NAME)

staticlib: $(LIBS)

libtinyxml2.a: tinyxml2.o
	$(AR) $(ARFLAGS)s $@ $^

clean :
	rm $(EXEC_NAME) $(OBJ_FILES) 

$(EXEC_NAME) : $(OBJ_FILES)
	$(CC) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

install :
	make
