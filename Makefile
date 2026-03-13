# BSGS Turbo - Makefile
# High-performance BSGS key search tool

CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -flto -pthread -Wall -Wextra -Wno-unused-parameter
SRCDIR = src
BUILDDIR = build
TARGET = bsgs_turbo

SOURCES = $(SRCDIR)/main.cpp $(SRCDIR)/bsgs_engine.cpp
HEADERS = $(SRCDIR)/uint256.h $(SRCDIR)/secp256k1.h $(SRCDIR)/hashtable.h $(SRCDIR)/bsgs_engine.h

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) -o $(TARGET) $(SOURCES)

# Debug build
debug: CXXFLAGS = -std=c++17 -g -O0 -pthread -Wall -Wextra -Wno-unused-parameter -DDEBUG
debug: $(TARGET)

# Test build
test: $(SRCDIR)/test_math.cpp $(HEADERS)
	$(CXX) -std=c++17 -O2 -I$(SRCDIR) -o test_math $(SRCDIR)/test_math.cpp
	./test_math

clean:
	rm -f $(TARGET) test_math *.o

# Windows-specific target (MinGW)
windows: $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) -o $(TARGET).exe $(SOURCES) -static
