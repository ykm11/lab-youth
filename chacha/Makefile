TARGET=main
CXXFLAG=-O3 -Wall -Wextra
LDFLAGS=-lgmpxx -lgmp
INCLUDES= -I ./include

all: $(TARGET)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAG) $(INCLUDES)

main: main.o chacha_poly.o
	$(CXX) main.o chacha_poly.o -o $@ $(CXXFLAG) $(LDFLAGS) $(INCLUDES)

test: test.cpp chacha_poly.o
	$(CXX) test.cpp chacha_poly.o -o $@ $(CXXFLAG) $(LDFLAGS) $(INCLUDES)
	./test

benchmark: benchmark.cpp chacha_poly.o
	$(CXX) benchmark.cpp chacha_poly.o -o $@ $(CXXFLAG) $(LDFLAGS) $(INCLUDES)

clean:
	$(RM) *.o ./test ./benchmark ./main

.PHONY:
	clean
	benchmark
	test
