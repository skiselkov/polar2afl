ACFUTILS=../libacfutils

CFLAGS=-g -O0 -I$(ACFUTILS)/src \
	$(shell $(ACFUTILS)/pkg-config-deps linux-64 --cflags)

LIBS=-L$(ACFUTILS)/qmake/lin64 -lacfutils \
	$(shell $(ACFUTILS)/pkg-config-deps linux-64 --libs) -lm

all : polar2afl

clean :
	rm -f polar2afl

polar2afl : polar2afl.c
	$(CC) -W -Wall -Werror $(CFLAGS) -o $@ $^ $(LIBS)
