COMPILE_FLAGS := g++ -g -std=c++17 -I ./include

makebuild: cartsweeper

SOURCES := $(wildcard ./src/*.cpp)
OBJECTS	:= $(patsubst ./src/%.cpp, ./bin/%.o, $(SOURCES))

cartsweeper: $(OBJECTS)
	@echo Linking...
	g++ -o ./cartsweeper $^

./bin/%.o: ./src/%.cpp
	@echo Compiling $< ...
	@$(COMPILE_FLAGS) -c $< -o $@

clean:
	@echo Cleaning...
	@rm cartsweeper
	@rm ./bin/*.o
