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

.PHONY: 
	clean 
	test
	test256
	test521
	benchmark
	benchmark256
	benchmark521
	benchmark521_mi

.SECONDARY: $(DEP)


./ec_test: tests/test.cpp fp256.o curve.o
	$(CXX) -o ec_test tests/test.cpp fp256.o curve.o $(CXXFLAGS) $(LDFLAGS)  -I ./include

./ec_test_mi: tests/test.cpp fp256.o curve.o
	$(CXX) -o ec_test_mi tests/test.cpp fp256.o curve.o $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -DUSE_MIMALLOC

test: 
	make test256
	make test521
	@echo "--- Elliptic Curve 256 bits TEST ---"
	./ec_test256
	@echo "--- Elliptic Curve 521 bits TEST ---"
	./ec_test521

./ec_benchmark: benchmarks/benchmark.cpp fp256.o curve.o
	$(CXX) -o ec_benchmark $? $(CXXFLAGS) $(LDFLAGS) -I ./include

./ec_benchmark_mi: benchmarks/benchmark.cpp fp256.o curve.o
	$(CXX) -o ec_benchmark_mi $? $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -DUSE_MIMALLOC

benchmark: 
	make benchmark256
	make benchmark521

# このあたりも綺麗にしたい
test256:
	$(CXX) -o fp256.o -c src/fp256.cpp $(CXXFLAGS) -I ./include
	$(CXX) -o curve.o -c src/curve.cpp $(CXXFLAGS) -I ./include
	$(CXX) -o ec_test256 tests/test.cpp fp256.o curve.o $(CXXFLAGS) $(LDFLAGS) -I ./include

test521:
	$(CXX) -o fp521.o -c src/fp521.cpp $(CXXFLAGS) -DYKM_ECC_SECP521 -I ./include
	$(CXX) -o curve.o -c src/curve.cpp $(CXXFLAGS) -DYKM_ECC_SECP521 -I ./include
	$(CXX) -o ec_test521 tests/test_secp521.cpp fp521.o curve.o $(CXXFLAGS) $(LDFLAGS) -DSECP521 -I ./include

benchmark256:
	$(CXX) -o fp256.o -c src/fp256.cpp $(CXXFLAGS) -I ./include
	$(CXX) -o curve.o -c src/curve.cpp $(CXXFLAGS) -I ./include
	$(CXX) -o ec_benchmark256 benchmarks/benchmark.cpp fp256.o curve.o $(CXXFLAGS) $(LDFLAGS) -I ./include

benchmark521:
	$(CXX) -o fp521.o -c src/fp521.cpp $(CXXFLAGS) -DYKM_ECC_SECP521 -I ./include
	$(CXX) -o curve.o -c src/curve.cpp $(CXXFLAGS) -DYKM_ECC_SECP521 -I ./include
	$(CXX) -o ec_benchmark521 benchmarks/benchmarkSecp521.cpp fp521.o curve.o $(CXXFLAGS) $(LDFLAGS) -DYKM_ECC_SECP521 -I ./include

benchmark521_mi:
	$(CXX) -o fp521.o -c src/fp521.cpp $(CXXFLAGS) -I ./include
	$(CXX) -o curve.o -c src/curve.cpp $(CXXFLAGS)  -I ./include
	$(CXX) -o ec_benchmark521_mi benchmarks/benchmarkSecp521.cpp fp521.o curve.o -DUSE_MIMALLOC $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) 


testEd25519: tests/test_twistedEd.cpp fp256.o curve.o
	$(CXX) -o ec_ed_test tests/test_twistedEd.cpp fp256.o curve.o $(CXXFLAGS) $(LDFLAGS)  -I ./include

benchmarkEd25519:
	$(CXX) -o fp256.o -c src/fp256.cpp $(CXXFLAGS) -I ./include
	$(CXX) -o curve.o -c src/curve.cpp $(CXXFLAGS) -I ./include
	$(CXX) -o ec_benchmarkEd25519 benchmarks/bench_ed25519.cpp fp256.o curve.o $(CXXFLAGS) $(LDFLAGS) -I ./include

