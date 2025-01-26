# Compiler and flags
CXX = /opt/homebrew/opt/llvm/bin/clang++
CXXFLAGS = -fcolor-diagnostics -fansi-escape-codes -g -pedantic-errors -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion -fopenmp -std=c++23 -O0 -arch arm64
INCLUDES = -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -I/opt/homebrew/Cellar/sfml/3.0.0/include -I/opt/homebrew/Cellar/SDL2/2.30.11/include -I/opt/homebrew/Cellar/glew/2.2.0_1/include -I/opt/homebrew/Cellar/glfw/3.4/include
LDFLAGS = -L/opt/homebrew/Cellar/sfml/3.0.0/lib -lsfml-graphics -lsfml-window -lsfml-system -L/opt/homebrew/Cellar/SDL2/2.30.11/lib -L/opt/homebrew/Cellar/glew/2.2.0_1/lib -L/opt/homebrew/Cellar/glfw/3.4/lib
LINKER_FLAGS = -lglfw -framework OpenGL -lglew

# Directories
SRC_DIR   = src
BUILD_DIR = build
TARGET    = $(BUILD_DIR)/main

# Source and object files
SRCS = $(wildcard $(SRC_DIR)/**/*.cpp) $(wildcard main.cpp)
OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)

# Default target
all: $(TARGET)

# Target to build the project
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LINKER_FLAGS) $(LDFLAGS) -o $(TARGET) $(OBJS)

# Rule to compile .cpp files into .o files
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean the build directory
clean:
	rm -rf $(BUILD_DIR)

.PHONY: clean
