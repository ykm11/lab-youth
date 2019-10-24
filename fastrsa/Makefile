TARGET=fastrsa
CXXFLAGS=-Wall -Wextra -O3 -std=c++17 -DNDEBUG
LFLAGS=-lgmpxx -lgmp


all: $(TARGET)


$(TARGET): main.cpp
	$(CXX) -o $< $(CXXFLAGS) $(LFLAGS)

rsa.o: rsa.cpp rsa.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)


.PHONY: benchmark clean

clean:
	rm *.o

./benchmark: test.cpp rsa.o
	$(CXX) -o $@ test.cpp rsa.o $(CXXFLAGS) $(LFLAGS)
	./benchmark
