# Compile up to 8 things at once
MAKEFLAGS = -j 8

# Compiler options
# CPPFLAGS = -Wall -Werror -pedantic
CPPFLAGS += -g
CPPFLAGS += -O3
CPPFLAGS += -std=c++11
#CPPFLAGS += -Wno-long-long

# Enable code profiling
#CPPFLAGS += -pg

HAMERLY_DIR = hamerly
SRC_HAMERLY = $(wildcard $(HAMERLY_DIR)/*.cpp)
OBJ_HAMERLY = $(SRC_HAMERLY:.cpp=.o)

SRC_COMMOM = $(filter-out HGMeans.cpp HGWrapper.cpp, $(wildcard *.cpp))
OBJ_COMMOM = $(SRC_COMMOM:.cpp=.o)

SRC_HG = HGMeans.cpp
OBJ_HG = $(SRC_HG:.cpp=.o)

all: hgmeans

hgmeans: $(OBJ_HAMERLY) $(OBJ_COMMOM) $(OBJ_HG)
	g++ $(OBJ_HAMERLY) $(OBJ_COMMOM) $(OBJ_HG) $(CPPFLAGS) $(LDFLAGS) -o hgmeans

.PHONY: clean all profile

profile:
	gprof hgmeans | less

clean:
	rm -f hgmeans $(OBJ_HAMERLY) $(OBJ_COMMOM) $(OBJ_HG) gmon.out