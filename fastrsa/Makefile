SRC=rsa.cpp benchmark.cpp
OBJ=$(SRC:.cpp=.o)

TARGET=rsa_benchmark
CXXFLAGS=-Wall -Wextra -O3 -std=c++17 -DNDEBUG
LDFLAGS=-lgmpxx -lgmp


all: $(TARGET)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $? $(LDFLAGS)

.PHONY: clean

clean:
	rm *.o $(TARGET)

