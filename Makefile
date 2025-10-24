CC=g++ #mpic++
CFLAGS=-std=c++17 -I /usr/local/include/eigen
OBJDIR=objs
LIBDIR=lib
SRCDIR=src
BINDINGDIR=python_bindings
SOURCES=$(SRCDIR)/main.cpp 
LIBSOURCES=$(LIBDIR)/CppInterfaceClass.cpp $(LIBDIR)/SolverClass.cpp $(LIBDIR)/PSRClass.cpp $(LIBDIR)/FixedHeightModelClass.cpp $(LIBDIR)/constants.cpp $(LIBDIR)/functions.cpp $(LIBDIR)/integrator.cpp $(LIBDIR)/read_write.cpp
BINDINGSOURCES=$(BINDINGDIR)/interface_class_binding.cpp 
LIBOBJECTS=$(LIBSOURCES:$(LIBDIR)/%.cpp=$(OBJDIR)/%.o)
SRCOBJECTS=$(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
BINDINGOBJECTS=$(BINDINGSOURCES:$(BINDINGDIR)/%.cpp=$(OBJDIR)/%.o)

EXECUTABLE=bin/loki

all: $(SOURCES) $(EXECUTABLE)

# binding: bin/libmylib.so $(BINDINGSOURCES)
# 	$(CC) $(CFLAGS) -shared -fPIC $$(python3 -m pybind11 --includes) -L bin -l mylib $(BINDINGSOURCES) -o $(BINDINGDIR)/cpp_interface$$(python3-config --extension-suffix)

whateverittakes: $(BINDINGSOURCES) $(LIBSOURCES)
	$(CC) $(CFLAGS) -shared -fPIC $$(python3 -m pybind11 --includes) $(BINDINGSOURCES) $(LIBSOURCES) -o $(BINDINGDIR)/cpp_interface$$(python3-config --extension-suffix)

lib: $(LIBSOURCES)
	$(CC) $(CFLAGS) -shared -fPIC $(LIBSOURCES) -o bin/libmylib.so


$(EXECUTABLE): $(LIBOBJECTS) $(SRCOBJECTS) 
	$(CC) $(CFLAGS) $(LIBOBJECTS) $(SRCOBJECTS) -o $@

$(OBJDIR)/%.o: $(LIBDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(BINDINGDIR)/%.cpp
	$(CC) -shared -fPIC $$(python3 -m pybind11 --includes) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o
