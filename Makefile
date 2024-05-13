# The compiler to be used
CXX = g++-14

# The list of flags passed to the compiler
CXXFLAGS = -O3 -Wall -std=c++11 -I ~/toofus `pkg-config --cflags --libs libtiff-4` -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 -lpthread

GLFLAGS = `pkg-config --cflags --libs glut` -framework OpenGL 


.PHONY: all clean

all: peanut

clean:
	rm -f *.o peanut

peanut: peanut.cpp 
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) peanut.cpp -o peanut

# The application that visualizes the conf files
see: see.cpp
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -o $@ see.cpp -Wno-deprecated $(GLFLAGS)
