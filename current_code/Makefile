CXX      ?= g++
CXXFLAGS ?= -O3 -fopenmp -march=native -std=c++17 -Wall
LDFLAGS  ?= -fopenmp

.PHONY: all clean

all: GL-panel

GL-panel: GL-panel.cpp common.hpp
	$(CXX) $(CXXFLAGS) -o $@ GL-panel.cpp $(LDFLAGS)

clean:
	rm -f GL-panel
