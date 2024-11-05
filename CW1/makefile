# Compiler
CXX = g++
CXXFLAGS = -Wall -std=c++11

# Target executable name
TARGET = your_program

# Source files
SOURCES = main.cpp matrix.cpp ode_solver.cpp opt_alg.cpp solution.cpp user_funs.cpp

# Header files
HEADERS = matrix.h ode_solver.h opt_alg.h solution.h user_funs.h

# Object files (derived from sources)
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Compile each source file into an object file
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJECTS) $(TARGET)

# Run target
run: $(TARGET)
	./$(TARGET)

# Phony targets
.PHONY: all clean run

