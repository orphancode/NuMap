CC = g++

CCFLAGS = -Wall -ansi -pedantic -O3 -funroll-loops
OBJECTS := $(patsubst %.cpp, %.o,$(wildcard *.cpp))

COMMON_HEADERS = ./
INCLUDE := $(COMMON_HEADERS)

all: $(OBJECTS)

%.o : %.cpp
	$(CC) -I$(INCLUDE) $(CCFLAGS) -c $<

clean: 
	rm -rf *.o 
