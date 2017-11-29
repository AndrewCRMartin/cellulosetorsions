CC=/usr/bin/cc
COPTS=-ansi -pedantic -Wall -L$(HOME)/lib -I$(HOME)/include -g
OFILES=cellulosetorsions.o
EXE=cellulosetorsions


$(EXE) : $(OFILES)
	$(CC) $(COPTS) -o $@ $< -lbiop -lgen -lm -lxml2

.c.o :
	$(CC) $(COPTS) -c -o $@ $<

clean :
	\rm -f $(OFILES)

distclean : clean
	\rm -f $(EXE)
