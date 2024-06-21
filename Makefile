# Compiler & Linker settings
CXX = g++
CXXFLAGS = -I ./inc -std=c++11
OPTFLAGS = -march=native -flto -funroll-loops -finline-functions -ffast-math -O3
WARNINGS = -g -Wall

# Valgrind for memory issue
CHECKCC = valgrind
CHECKFLAGS = --leak-check=full -s --show-leak-kinds=all --track-origins=yes 

# Cppcheck for static analysis
CPPCHECKCC = cppcheck
CPPCHECKFLAGS = --enable=all --inconclusive --std=c++11 --language=c++ --suppress=missingIncludeSystem --suppress=style --suppress=unusedFunction --suppress=unreadVariable --suppress=unmatchedSuppression --suppress=unusedStructMember --suppress=functionConst

# Source files and object files
SRCDIR = src
OBJDIR = obj
INCDIR = inc
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
DEPS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.d,$(SRCS))

# Control the build verbosity
ifeq ("$(VERBOSE)","1")
    Q :=
    VECHO = @true
else
    Q := @
    VECHO = @printf
endif

.PHONY: all check clean calc

# Name of the executable
BIN = cadb0015

all: $(BIN)

$(OBJDIR):
	@mkdir $(OBJDIR)

$(BIN): main.cpp $(OBJS)
	$(VECHO) "	LD\t$@\n"
	$(Q)$(CXX) $(WARNINGS) $(CXXFLAGS) $(OPTFLAGS) $^ -o $@ $(LINKER)

# Include generated dependency files
-include $(DEPS)

# Compilation rule for object files with automatic dependency generation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR) Makefile
	$(VECHO) "	CC\t$@\n"
	$(Q)$(CXX) $(WARNINGS) $(CXXFLAGS) $(OPTFLAGS) -MMD -c $< -o $@

run1:
	./$(BIN) testcase/sampleCase

run2:
	./$(BIN) testcase/sample.txt

run3:
	./$(BIN) testcase/testcase1.txt

check:
	$(CHECKCC) $(CHECKFLAGS) ./$(BIN) testcase/testcase1.txt

cppcheck:
	$(CPPCHECKCC) $(CPPCHECKFLAGS) -I$(INCDIR) $(SRCDIR) main.cpp $(INCDIR)/*.h


clean:
	rm -rf $(OBJDIR) $(BIN)
