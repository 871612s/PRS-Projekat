# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++23 -O3 -Wall -Wextra

# Executable name
TARGET = projekat

# Source files
SRC = projekat.cpp

# Build target
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

# Remove generated files
clean:
	rm -f $(TARGET)