SRC=mr.cpp mr_benchmark.cpp
OBJ=$(SRC:.cpp=.o)

TARGET=mr_benchmark
LDFLAGS=-lgmpxx -lgmp
CXXFLAGS=-O3 -Wall -Wextra -std=c++17


all: $(TARGET)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) 

$(TARGET): $(OBJ) 
	$(CXX) -o $@ $? $(LDFLAGS)


.PHONY: clean benchmark

clean:
	$(RM) -rf $(OBJ) $(TARGET)

benchmark:
	make clean
	make 
	./mr_benchmark
