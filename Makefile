CC = gcc
CFLAGS = -fopenmp
TARGET = heat_serial
SRC = heat_serial.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)
