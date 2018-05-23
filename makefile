SRC_DIR = src
OBJ_DIR = $(SRC_DIR)/obj
TEST_DIR = $(SRC_DIR)/test
OBJ = $(addprefix $(OBJ_DIR)/, oisdifference.o fitshelper.o)
HEADERS = $(addprefix $(SRC_DIR)/, oisdifference.h fitshelper.h)
LIBS = -lm -lcfitsio
OPTS = -std=gnu99

all: oisdifference

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

oisdifference: $(SRC_DIR)/main.c $(OBJ) $(HEADERS)
	$(CC) $(OPTS) $(SRC_DIR)/main.c $(OBJ) -o oisdifference $(LIBS)

$(OBJ_DIR)/oisdifference.o: $(OBJ_DIR) $(SRC_DIR)/oisdifference.c $(HEADERS)
	$(CC) $(OPTS) -c $(SRC_DIR)/oisdifference.c -o $(OBJ_DIR)/oisdifference.o

$(OBJ_DIR)/fitshelper.o: $(OBJ_DIR) $(SRC_DIR)/fitshelper.c $(SRC_DIR)/fitshelper.h
	$(CC) $(OPTS) -c $(SRC_DIR)/fitshelper.c -o $(OBJ_DIR)/fitshelper.o

.PHONY: all clean install

install:
	mv oisdifference /usr/local/bin

test:
	$(CC) $(OPTS) $(OBJ) src/test/test_oisdifference.c -I$(SRC_DIR) -o $(TEST_DIR)/test -lm -lcfitsio
	cd $(TEST_DIR); ./test

clean:
	rm -f oisdifference
	rm -rf $(OBJ_DIR)
	rm -f $(TEST_DIR)/test
