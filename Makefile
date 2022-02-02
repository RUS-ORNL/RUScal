PROG = ruscal
SRCS = rus.c matrix.c
HDRS = matrix.h
OBJS = $(SRCS:.c=.o)
C = gcc
CFLAGS = -Wall -g 
LDFLAGS = -lm

$(PROG) : $(OBJS) $(HDRS)
	$(C) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(PROG)

clean: 
	rm  -f *~ $(PROG) $(OBJS)
