CC=cc
COPTS=-L$(HOME)/lib -I$(HOME)/include -O3
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
