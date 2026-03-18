CXX      := g++
CXXFLAGS := --std=c++23 -Wall -Wextra -fopenmp

TARGET   := probeFuzz
SRC      := probeFuzz.cpp

# Release (default)
RELEASE_FLAGS := -O3 -march=native

.PHONY: all
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(RELEASE_FLAGS) $(LDFLAGS) -o $@ $<

# Debug
DEBUG_FLAGS := -O0 -g3 -fsanitize=address,undefined

.PHONY: debug
debug: $(SRC)
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) $(LDFLAGS) -o $(TARGET)_debug $<

# Clean
.PHONY: clean
clean:
	rm -f $(TARGET) $(TARGET)_debug
