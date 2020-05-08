SRC=fp256.cpp curve.cpp main.cpp
OBJ=$(SRC:.cpp=.o)
DEP=$(SRC:.cpp=.d)

TARGET=main
CXXFLAGS=-Wall -Wextra -O3 -std=c++17 -DNDEBUG -DYKM_ECC_USE_MPN
LDFLAGS=-lgmpxx -lgmp -lpthread
INCLUDE=-I /usr/local/lib/mimalloc-1.6/include /usr/local/lib/mimalloc-1.6/libmimalloc.a -I ./include

all: $(TARGET)

%.o: src/%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) -MMD -MP -MF $(@:.o=.d) -I ./include

-include $(DEP)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $? $(LDFLAGS)

clean:
	$(RM) -rf $(OBJ) $(DEP) $(TARGET) ec_* fp*.o fp*.d

.PHONY: clean benchmark

.SECONDARY: $(DEP)


./ec_test: tests/test.cpp fp256.o curve.o
	$(CXX) -o ec_test tests/test.cpp fp256.o curve.o $(CXXFLAGS) $(LDFLAGS)  

./ec_test_mi: tests/test.cpp fp256.o curve.o
	$(CXX) -o ec_test_mi tests/test.cpp fp256.o curve.o $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -DUSE_MIMALLOC

test: ./ec_test_mi
	./ec_test_mi

./ec_benchmark: benchmarks/benchmark.cpp fp256.o curve.o
	$(CXX) -o ec_benchmark $? $(CXXFLAGS) $(LDFLAGS)  

./ec_benchmark_mi: benchmarks/benchmark.cpp fp256.o curve.o
	$(CXX) -o ec_benchmark_mi $? $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -DUSE_MIMALLOC

benchmark: ./ec_benchmark_mi
	./ec_benchmark_mi

# このあたりも綺麗にしたい
test521:
	$(CXX) -o fp521.o -c src/fp521.cpp $(CXXFLAGS) -DYKM_ECC_SECP521 -I ./include
	$(CXX) -o curve.o -c src/curve.cpp $(CXXFLAGS) -DYKM_ECC_SECP521 -I ./include
	$(CXX) -o ec_test521 tests/test_secp521.cpp fp521.o curve.o $(CXXFLAGS) $(LDFLAGS) -DSECP521 -I ./include

benchmark521:
	$(CXX) -o fp521.o -c src/fp521.cpp $(CXXFLAGS) -DYKM_ECC_SECP521 -I ./include
	$(CXX) -o curve.o -c src/curve.cpp $(CXXFLAGS) -DYKM_ECC_SECP521 -I ./include
	$(CXX) -o ec_benchmark521 benchmarks/benchmarkSecp521.cpp fp521.o curve.o $(CXXFLAGS) $(LDFLAGS) -DYKM_ECC_SECP521 $(INCLUDE)

benchmark521_mi:
	$(CXX) -o FP.o -c src/fp521.cpp $(CXXFLAGS) -I ./include
	$(CXX) -o curve.o -c src/curve.cpp $(CXXFLAGS)  -I ./include
	$(CXX) -o ec_benchmark521_mi benchmarks/benchmarkSecp521.cpp FP.o curve.o $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) 

