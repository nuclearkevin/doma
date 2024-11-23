COMPILE_FLAGS := g++ -g -std=c++17 -I ./include

makebuild: doma

SOURCES := $(wildcard ./src/*.cpp)
OBJECTS	:= $(patsubst ./src/%.cpp, ./bin/%.o, $(SOURCES))

doma: $(OBJECTS)
	@echo Linking...
	g++ -o ./doma $^

./bin/%.o: ./src/%.cpp
	@echo Compiling $< ...
	@$(COMPILE_FLAGS) -c $< -o $@

clean:
	@echo Cleaning...
	@rm doma
	@rm ./bin/*.o
