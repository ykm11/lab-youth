SRC=FP.cpp curve.cpp main.cpp
OBJ=$(SRC:.cpp=.o)
DEP=$(SRC:.cpp=.d)

TARGET=main
CXXFLAGS=-Wall -Wextra -O3 -std=c++17 -DNDEBUG -DUSE_MPN
LDFLAGS=-lgmpxx -lgmp -lpthread
INCLUDE=-I /usr/local/lib/mimalloc-1.6/include /usr/local/lib/mimalloc-1.6/libmimalloc.a -I ./include

all: $(TARGET)

%.o: src/%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) -MMD -MP -MF $(@:.o=.d) -I ./include

-include $(DEP)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $? $(LDFLAGS)

clean:
	$(RM) -rf $(OBJ) $(DEP) $(TARGET) ec_*

.PHONY: clean benchmark

.SECONDARY: $(DEP)


./ec_test: test.cpp FP.o curve.o
	$(CXX) -o ec_test test.cpp FP.o curve.o $(CXXFLAGS) $(LDFLAGS)  

./ec_test_mi: test.cpp FP.o curve.o
	$(CXX) -o ec_test_mi test.cpp FP.o curve.o $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -DUSE_MIMALLOC

test: ./ec_test_mi
	./ec_test_mi

./ec_benchmark: benchmark.cpp FP.o curve.o
	$(CXX) -o ec_benchmark $? $(CXXFLAGS) $(LDFLAGS)  

./ec_benchmark_mi: benchmark.cpp FP.o curve.o
	$(CXX) -o ec_benchmark_mi $? $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -DUSE_MIMALLOC

benchmark: ./ec_benchmark_mi
	./ec_benchmark_mi

test521:
	$(CXX) -o FP.o -c FP.cpp $(CXXFLAGS) -DSECP521
	$(CXX) -o curve.o -c curve.cpp $(CXXFLAGS) -DSECP521
	$(CXX) -o ec_test521 test_secp521.cpp FP.o curve.o $(CXXFLAGS) $(LDFLAGS) -DSECP521

benchmark521:
	$(CXX) -o FP.o -c FP.cpp $(CXXFLAGS) -DSECP521
	$(CXX) -o curve.o -c curve.cpp $(CXXFLAGS) -DSECP521
	$(CXX) -o ec_benchmark521 benchmarkSecp521.cpp FP.o curve.o $(CXXFLAGS) $(LDFLAGS) -DSECP521

benchmark521_mi:
	$(CXX) -o FP.o -c FP.cpp $(CXXFLAGS) 
	$(CXX) -o curve.o -c curve.cpp $(CXXFLAGS) 
	$(CXX) -o ec_benchmark521_mi benchmarkSecp521.cpp FP.o curve.o $(CXXFLAGS) $(LDFLAGS) $(INCLUDE)

