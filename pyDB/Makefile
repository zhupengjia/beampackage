CC=gcc
SFLAGS=-fPIC -shared
SRC=$(wildcard bpm*.c)
OBJS=$(SRC:.c=.so)

bpm:$(OBJS)

%.so:%.c
	$(CC) $< -o $@ $(SFLAGS)

clean:
	rm *.so -f
