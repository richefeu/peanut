# The compiler to be used
CXX = g++-13

# The list of flags passed to the compiler
CXXFLAGS = -O3 -Wall -std=c++11 -I ~/toofus `pkg-config --cflags --libs libtiff-4` -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 -lpthread  

.PHONY: all clean

all: peanut

clean:
	rm -f *.o peanut

peanut: peanut.cpp 
	$(CXX) $(CXXFLAGS) peanut.cpp -o peanut
