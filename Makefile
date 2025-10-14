# FILES = src/*.cpp lib/*.cpp

# CC = mpic++
# CFLAGS = -std=c++17 -I /usr/local/include/eigen

# all:
# 	mkdir -p bin
# 	${CC} ${CFLAGS} ${FILES} -o bin/loki

# clean:
# 	rm bin/loki


CC=g++ #mpic++
CFLAGS=-std=c++17 -I /usr/local/include/eigen
OBJDIR=objs
LIBDIR=lib
SRCDIR=src
SOURCES=$(SRCDIR)/main_test.cpp 
LIBSOURCES=$(LIBDIR)/CppInterfaceClass.cpp $(LIBDIR)/SolverClass.cpp $(LIBDIR)/PSRClass.cpp $(LIBDIR)/FixedHeightModelClass.cpp $(LIBDIR)/constants.cpp $(LIBDIR)/functions.cpp $(LIBDIR)/integrator.cpp
OBJECTS=$(LIBSOURCES:$(LIBDIR)/%.cpp=$(OBJDIR)/%.o)
OBJECTS+=$(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
EXECUTABLE=bin/test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) -o $@

$(OBJDIR)/%.o: $(LIBDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o
