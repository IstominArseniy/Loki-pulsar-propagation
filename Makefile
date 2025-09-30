FILES = src/*.cpp lib/*.cpp

CC = mpic++
CFLAGS = -std=c++17 -I /usr/local/include/eigen

all:
	mkdir -p bin
	${CC} ${CFLAGS} ${FILES} -o bin/loki

clean:
	rm bin/loki
