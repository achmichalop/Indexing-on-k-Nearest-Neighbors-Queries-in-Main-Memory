OS := $(shell uname)
ifeq ($(OS),Darwin)
		CC      = /usr/local/bin/g++-10
        CFLAGS  = -O3 -g -mavx -std=c++14 -w -march=native -I/usr/local/opt/libomp/include -I/usr/local/include/ -fopenmp 
        LDFLAGS = -L/usr/local/opt/libomp/lib   -L/usr/local/lib/  
else
        CC      = g++
        CFLAGS  = -O3 -g -mavx -std=c++14 -w -march=native -fopenmp -I/usr/local/include/ 
        LDFLAGS = -L/usr/local/lib/
endif


SOURCES = containers/relation.cpp PriorityQueue/priorityQueue.cpp quering/process.cpp
OBJECTS = $(SOURCES:.cpp=.o)
	
	
all: twoLevel 

twoLevel: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) main_twolevel.cpp -o twoLevel $(LDADD)

	
.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

.cc.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf containers/*.o
	rm -rf partitioning/*.o
	rm -f twoLevel
	rm -f PriorityQueue/*.o
	rm -f quering/*.o
