CXX=g++         # C++ compiler
CXXFLAGS=-Wall  # compiler flags

# source files
SRC=settings.cpp main.cpp


all: debug release

debug:
	$(CXX) $(SRC) $(CXXFLAGS) -g -o test2.obj

release:
	$(CXX) $(SRC) $(CXXFLAGS) -Ofast -DNDEBUG -o test2_release

clean:
	rm test2 test2_release