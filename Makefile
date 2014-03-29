CC=g++
SFLAGS=-fPIC -shared
SRC=$(wildcard *.c)
OBJS=$(SRC:.c=.so)

all:$(OBJS)

%.so:%.c
	$(CC) $< -o $@ $(SFLAGS)

clean:
	rm *.so -f
