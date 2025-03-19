FILES = src/*.cpp lib/*.cpp

CC = g++
CFLAGS = -std=c++14

all:
	mkdir -p bin
	${CC} ${CFLAGS} ${FILES} -I /usr/local/include/eigen -o bin/loki

clean:
	rm bin/loki
