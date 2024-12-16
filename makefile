COMPILE_FLAGS := g++ -O2 -std=c++17 -I ./include -I ./contrib/pugixml/src -I ./contrib/argparse/include

makebuild: dir pugi doma

SOURCES := $(wildcard ./src/*.cpp)
OBJECTS	:= $(patsubst ./src/%.cpp, ./bin/%.o, $(SOURCES))

PUGI_SOURCES := $(wildcard ./contrib/pugixml/src/*.cpp)
PUGI_OBJECTS := $(patsubst ./contrib/pugixml/src/%.cpp, ./bin/%.o, $(PUGI_SOURCES))

dir:
	@mkdir -p bin

pugi: ${PUGI_OBJECTS}
	@echo Building pugixml
	@$(COMPILE_FLAGS) -c $< -o $@

doma: $(OBJECTS)
	@echo Linking...
	g++ -o ./doma ${PUGI_OBJECTS} $^

./bin/%.o: ./src/%.cpp
	@echo Compiling $< ...
	@$(COMPILE_FLAGS) -c $< -o $@

./bin/%.o: ./contrib/pugixml/src/%.cpp
	@echo Compiling $< ...
	@$(COMPILE_FLAGS) -c $< -o $@

clean:
	@echo Cleaning...
	@rm ./bin/*.o
	@rm doma
