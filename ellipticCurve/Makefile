SRC=FP.cpp curve.cpp main.cpp
OBJ=$(SRC:.cpp=.o)
DEP=$(SRC:.cpp=.d)

TARGET=main
CXXFLAGS=-Wall -Wextra -O3 -std=c++17 -DNDEBUG -DUSE_MPN
LDFLAGS=-lgmpxx -lgmp -lpthread

all: $(TARGET)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) -MMD -MP -MF $(@:.o=.d)

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
	$(CXX) -o ec_test_mi test.cpp FP.o curve.o $(CXXFLAGS) $(LDFLAGS) -I /usr/local/lib/mimalloc-1.6/include /usr/local/lib/mimalloc-1.6/libmimalloc.a -DUSE_MIMALLOC

test: ./ec_test_mi
	./ec_test_mi


./ec_benchmark: test.cpp FP.o curve.o
	$(CXX) -o ec_benchmark benchmark.cpp FP.o curve.o $(CXXFLAGS) $(LDFLAGS)  

./ec_benchmark_mi: test.cpp FP.o curve.o
	$(CXX) -o ec_benchmark_mi benchmark.cpp FP.o curve.o $(CXXFLAGS) $(LDFLAGS) -I /usr/local/lib/mimalloc-1.6/include /usr/local/lib/mimalloc-1.6/libmimalloc.a -DUSE_MIMALLOC

benchmark: ./ec_benchmark_mi
	./ec_benchmark_mi

