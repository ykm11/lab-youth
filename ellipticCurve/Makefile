target=main
CXXFLAGS=-Wall -Wextra -O3 -std=c++17
LDFLAGS=-lgmpxx -lgmp

./main: main.cpp FP.h curve.h
	$(CXX) -o main main.cpp -lgmpxx -lgmp -std=c++17 

clean:
	rm ./ec_test
	rm ./main

./ec_test: test.cpp FP.h
	$(CXX) -o ec_test test.cpp $(CXXFLAGS) $(LDFLAGS)  

test: ./ec_test
	./ec_test
