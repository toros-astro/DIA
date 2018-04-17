CC = gcc
SRC_DIR = src
OBJ_DIR = $(SRC_DIR)/obj
OBJ = $(addprefix $(OBJ_DIR)/, oisdifference.o fitshelper.o)
HEADERS = $(addprefix $(SRC_DIR)/, oisdifference.h fitshelper.h)
LIBS = -lm -lcfitsio

all: oisdifference

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

oisdifference: $(SRC_DIR)/main.c $(OBJ) $(HEADERS)
	$(CC) $(SRC_DIR)/main.c $(OBJ) -o oisdifference $(LIBS)

$(OBJ_DIR)/oisdifference.o: $(OBJ_DIR) $(SRC_DIR)/oisdifference.c $(HEADERS)
	$(CC) -c $(SRC_DIR)/oisdifference.c -o $(OBJ_DIR)/oisdifference.o

$(OBJ_DIR)/fitshelper.o: $(OBJ_DIR) $(SRC_DIR)/fitshelper.c $(SRC_DIR)/fitshelper.h
	$(CC) -c $(SRC_DIR)/fitshelper.c -o $(OBJ_DIR)/fitshelper.o

.PHONY: all clean install

install:
	mv oisdifference /usr/local/bin

clean:
	rm -f oisdifference
	rm -rf $(OBJ_DIR)
