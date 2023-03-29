# Makefile for building spkmeans executable

# Compiler options
CC = gcc
CFLAGS = -Wall -g -Wextra -Werror -pedantic-errors
LDFLAGS = -lm 
# Source files
SOURCES = spkmeans.c main.c
HEADERS = spkmeans.h

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Build target
TARGET = spkmeans

# Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@ 

clean:
	rm -f $(TARGET) $(OBJECTS)