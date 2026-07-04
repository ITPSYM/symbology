# compiler
CXX := g++
CXXFLAGS := -O3 -flto -march=native -mtune=native -std=c++20 -I.
LDLIBS := -lflint -lgmp -lmimalloc

ifeq ($(OS),Windows_NT)
EXE := bootstrap.exe
LDLIBS += -ltbb12
else
EXE := bootstrap
LDLIBS += -ltbb
endif

# target
all: $(EXE)

# executable file
$(EXE): bootstrap.cpp bootstrap.hpp projection.hpp
	$(CXX) bootstrap.cpp -o $@ $(CXXFLAGS) $(LDLIBS)

# clean target
clean:
	rm -f bootstrap bootstrap.exe
