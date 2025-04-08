CC = gcc
CFLAGS = -fopenmp
TARGET_SERIAL = heat_serial
TARGET_PARALLEL = heat_parallel
SRC_SERIAL = heat_serial.c
SRC_PARALLEL = heat_parallel.c

all: $(TARGET_SERIAL) $(TARGET_PARALLEL)

serial: $(TARGET_SERIAL)

parallel: $(TARGET_PARALLEL)

$(TARGET_SERIAL): $(SRC_SERIAL)
	$(CC) -o $(TARGET_SERIAL) $(SRC_SERIAL)

$(TARGET_PARALLEL): $(SRC_PARALLEL)
	$(CC) $(CFLAGS) -o $(TARGET_PARALLEL) $(SRC_PARALLEL)

clean:
	rm -f $(TARGET_SERIAL) $(TARGET_PARALLEL)
	rm -f *.bmp
