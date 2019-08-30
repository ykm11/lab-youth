SRC=main.cpp
OBJ=$(SRC:.cpp=.o)
DEP=$(SRC:.cpp=.d)

TARGET=main
CXXFLAGS=-Wall -Wextra -O3 -std=c++17
LDFLAGS=-lgmpxx -lgmp

all: $(TARGET)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) -MMD -MP -MF $(@:.o=.d)

-include $(DEP)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $? $(LDFLAGS)

clean:
	$(RM) -rf $(OBJ) $(DEP) $(TARGET)

.PHONY: clean

.SECONDARY: $(DEP)


./ec_test: test.cpp FP.h curve.h
	$(CXX) -o ec_test test.cpp $(CXXFLAGS) $(LDFLAGS)  

test: ./ec_test
	./ec_test
