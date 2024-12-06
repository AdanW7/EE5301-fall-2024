# C++ compilation flags
CFLAGS=-O3
CPP = g++

# Python virtual environment setup
VENV = .venv
PYTHON = $(VENV)/bin/python
PIP = $(PYTHON) -m pip
REQUIREMENTS = requirements.txt

# Targets for both C++ and Python
.PHONY: all clean install run

# Default target
all: placer main install

# C++ targets
placer: suraj_parser.o main.o
	$(CPP) $(CFLAGS) -o placer suraj_parser.o main.o

main: suraj_parser.o main.o
	$(CPP) $(CFLAGS) -o main suraj_parser.o main.o

suraj_parser.o: suraj_parser.cpp
	$(CPP) $(CFLAGS) -c suraj_parser.cpp

main.o: main.cpp suraj_parser.h
	$(CPP) $(CFLAGS) -c main.cpp

# Python virtual environment and requirements installation
$(VENV)/bin/activate: 
	python3 -m venv $(VENV)

install: $(VENV)/bin/activate $(REQUIREMENTS)
	$(PIP) install --upgrade pip
	$(PIP) install -r $(REQUIREMENTS)

# Clean up build files and virtual environment
clean: 
	rm -rf *.o *.csv *.png placer main requirements.txt
	rm -rf $(VENV)

# Run Python script (graph_maker.py)
#	make run file=toy01 
#	make run file=toy02 
#	make run file=ibm01 
run:
	$(PYTHON) graph_maker.py $(file)

# If no requirements.txt exists, create it
requirements.txt:
	echo "matplotlib" > $(REQUIREMENTS)
	echo "pandas" >> $(REQUIREMENTS)
