target=main
CXXFLAGS=-Wall -Wextra -O3 -std=c++17
LDFLAGS=-lgmpxx -lgmp

all:
	g++ ${target}.cpp -lgmpxx -lgmp -o $(target) -std=c++17 

clean:
	-rm test

./ec_test: 
	$(CXX) -o ec_test test.cpp $(CXXFLAGS) $(LDFLAGS)  

test: ./ec_test
	./ec_test
