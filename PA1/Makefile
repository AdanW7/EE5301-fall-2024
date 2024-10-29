CC=g++
TARGET=sta
CFLAGS=-Wall -std=c++11 -g


SRCS = main.cpp Circuit.cpp CircuitNode.cpp GateDatabase.cpp
OBJS = $(SRCS:.cpp=.o)

# Default target to compile the program
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

# Compile .cpp files into .o object files
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Run the executable
run: $(TARGET)
	./$(TARGET)

# Clean the directory by removing object files and the executable
clean:
	@echo "Cleaning"
	rm -f $(OBJS) $(TARGET)

# Test target
test: $(TARGET)
	@echo "Testing"
	./sta test/NLDM_lib_max2Inp test/cleaned_iscas89_99_circuits/c17.isc

.PHONY: all clean test run
