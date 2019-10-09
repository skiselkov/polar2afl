# CDDL HEADER START
#
# This file and its contents are supplied under the terms of the
# Common Development and Distribution License ("CDDL"), version 1.0.
# You may only use this file in accordance with the terms of version
# 1.0 of the CDDL.
#
# A full copy of the text of the CDDL should have accompanied this
# source.  A copy of the CDDL is also available via the Internet at
# http://www.illumos.org/license/CDDL.
#
# CDDL HEADER END
#
# Copyright 2019 Saso Kiselkov. All rights reserved.

OBJS = polar2afl.o avl.o
OUTBIN=polar2afl

CFLAGS=-W -Wall -Werror -g -O2
LIBS=-lm

ifeq ($(findstring mingw,$(CC)),mingw)
	OUTBIN := $(OUTBIN).exe
endif

all : $(OUTBIN)

clean :
	rm -f $(OUTBIN) $(OBJS)

$(OUTBIN) : $(OBJS)
	$(CC) -W -Wall -Werror -o $@ $^ $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $^
