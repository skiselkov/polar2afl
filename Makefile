OBJS = polar2afl.o avl.o

CFLAGS=-W -Wall -Werror -g -O2
LIBS=-lm

all : polar2afl

clean :
	rm -f polar2afl $(OBJS)

polar2afl : $(OBJS)
	$(CC) -W -Wall -Werror -o $@ $^ $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $^
