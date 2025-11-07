CC=g++ #mpic++
CCMPI=mpic++
CFLAGS=-std=c++17 -I /usr/local/include/eigen
OBJDIR=bin/objs
LIBDIR=lib
SRCDIR=src
BINDINGDIR=python_bindings
SOURCES=$(SRCDIR)/main.cpp 
SOURCES_MPI=$(SRCDIR)/main_mpi.cpp 
LIBSOURCES=$(LIBDIR)/CppInterfaceClass.cpp $(LIBDIR)/SolverClass.cpp $(LIBDIR)/PSRClass.cpp $(LIBDIR)/FixedHeightModelClass.cpp $(LIBDIR)/constants.cpp $(LIBDIR)/functions.cpp $(LIBDIR)/integrator.cpp $(LIBDIR)/read_write.cpp
BINDINGSOURCES=$(BINDINGDIR)/interface_class_binding.cpp 
LIBOBJECTS=$(LIBSOURCES:$(LIBDIR)/%.cpp=$(OBJDIR)/%.o)
SRCOBJECTS=$(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
BINDINGOBJECTS=$(BINDINGSOURCES:$(BINDINGDIR)/%.cpp=$(OBJDIR)/%.o)

EXECUTABLE=bin/loki
EXECUTABLE_MPI=bin/loki_mpi


all: main main_mpi

main: $(SOURCES) $(EXECUTABLE)

main_mpi: $(SOURCES_MPI) $(EXECUTABLE_MPI)

whateverittakes: $(BINDINGSOURCES) $(LIBSOURCES)
	$(CC) $(CFLAGS) -shared -fPIC $$(python3 -m pybind11 --includes) $(BINDINGSOURCES) $(LIBSOURCES) -o bin/cpp_interface$$(python3-config --extension-suffix)

lib: $(LIBSOURCES)
	$(CC) $(CFLAGS) -shared -fPIC $(LIBSOURCES) -o bin/libmylib.so

$(EXECUTABLE): $(LIBOBJECTS) $(OBJDIR)/main.o 
	$(CC) $(CFLAGS) $(LIBOBJECTS) $(OBJDIR)/main.o -o $@

$(EXECUTABLE_MPI): $(LIBOBJECTS) $(OBJDIR)/main_mpi.o 
	$(CCMPI) $(CFLAGS) $(LIBOBJECTS) $(OBJDIR)/main_mpi.o -o $@	

$(OBJDIR)/%.o: $(LIBDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/main.o: $(SRCDIR)/main.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/main_mpi.o: $(SRCDIR)/main_mpi.cpp
	$(CCMPI) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(BINDINGDIR)/%.cpp
	$(CC) -shared -fPIC $$(python3 -m pybind11 --includes) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o
