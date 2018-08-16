OBJ_PATH = objects
SRC_PATH = source

SRCS = $(wildcard $(SRC_PATH)/*.cpp)

OBJS = $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(SRCS)))

CC = g++ -std=c++17 -O3 -g

all: jDE-demo

jDE-demo: $(OBJS)
	$(CC) $^ -o $@

$(OBJ_PATH)/%.o : $(SRC_PATH)/%.cpp
	$(CC) -o $@ -c $<

clean:
	-rm -f $(OBJ_PATH)/*.o jDE-demo

run:
	./jDE-demo
