TARGET=fastrsa
CXXFLAGS=-Wall -Wextra -O3 -std=c++17 -DNDEBUG
LFLAGS=-lgmpxx -lgmp


all: $(TARGET)


$(TARGET): main.cpp
	$(CXX) -o $< $(CXXFLAGS) $(LFLAGS)

.PHONY: benchmark


benchmark: test.cpp rsa.h
	$(CXX) -o ./benchmark test.cpp $(CXXFLAGS) $(LFLAGS)
	./benchmark
