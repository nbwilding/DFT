# MAKEFILE for C-program film
# OPTIONS

NAME     = DFT
INCLUDES = -I.
OUT      = \".\"
IN       = \".\"


CFLAGS  = -O3 $(INCLUDES) -DOUTDIR=$(OUT) -DINDIR=$(IN)  #-DNONAME
CC      =  gcc
LDFLAGS = -lm 

# suffix rule: .c needed for .o

.c.o  :  
	$(CC) $(CFLAGS) -c $*.c

#macro definition

DFT =	dft.o 

# The following commands (second line) produce the files (targets)
# listed on the left of the colon given the prerequisites on the RHS of the first line

DFT:		$(DFT)
		$(CC) -o $(NAME) $(CFLAGS) $(DFT) $(LDFLAGS)
#	        rm DFT.o

debug:		DFT.c makefile
		$(CC) $(CFLAGS) -g -c DFT.c

DFT.o:		DFT.c convlv.c numrec.c makefile
		$(CC) $(CFLAGS) -c DFT.c

new:
		touch *.c

#type "make clean" to execute this command

clean:
		rm *.o

