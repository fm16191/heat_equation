CCX=g++
CXXFLAGS=-g3 -Wall -Wextra

all: heat_1d

heat_1d:
	$(CCX) $(CXXFLAGS) heat_equation_1d.cxx -o $@
