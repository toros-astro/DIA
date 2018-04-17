CC = gcc
SRC_DIR = routines/C
LIBS = -lm -lcfitsio

all: oisdifference

oisdifference: $(SRC_DIR)/oisdifference.c
	$(CC) $(SRC_DIR)/oisdifference.c -o oisdifference $(LIBS)

.PHONY: all clean install

install:
	mv oisdifference /usr/local/bin

clean:
	rm oisdifference
