TARGET=main
OBJ=chacha_poly.o
CXXFLAG=-O3 -Wall -Wextra
LDFLAGS=-lgmpxx -lgmp
INCLUDES= -I ./include
SRCDIR=src

all: $(TARGET)

%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAG) $(INCLUDES)

main: main.o $(OBJ)
	$(CXX) main.o $(OBJ) -o $@ $(CXXFLAG) $(LDFLAGS) $(INCLUDES)

test: test.cpp $(OBJ)
	$(CXX) test.cpp $(OBJ) -o $@ $(CXXFLAG) $(LDFLAGS) $(INCLUDES)
	./test

benchmark: benchmark.cpp $(OBJ)
	$(CXX) benchmark.cpp $(OBJ) -o $@ $(CXXFLAG) $(LDFLAGS) $(INCLUDES)

clean:
	$(RM) *.o ./test ./benchmark ./main

.PHONY:
	clean
	benchmark
	test
