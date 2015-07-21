TARGET = vetospectrumroofit 
CC = g++
ROOTCINT = $(ROOTSYS)/bin/rootcint
SRCS = $(addsuffix .C, $(TARGET))
OBJS = $(addsuffix .o, $(notdir $(basename $(SRCS))))
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags) #-Wall -fPIC -g 
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs) -lProof -lProofPlayer -lRooFit -lRooFitCore
ROOTLDFLAGS = $(shell $(ROOTSYS)/bin/root-config --ldflags)


vetospectrumroofit: $(SRCS)
	$(CC) $(ROOTCFLAGS) $(ROOTLIBS) $^ -g -o $@ 

.PHONY: clean


clean:
	rm -rf $(TARGET)

