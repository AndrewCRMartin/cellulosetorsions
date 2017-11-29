CC     = cc
#CC     = x86_64-w64-mingw32-gcc  # For windows cross compiler

XML    = # -DXML_SUPPORT $(shell xml2-config --cflags)
XMLLIB = # -lxml2
#
GUNZIP = -DGUNZIP_SUPPORT
#
COPT   = -O3 $(XML) $(GUNZIP)
LIBS   = $(XMLLIB) -lm
LFILES = bioplib/BuildConect.o     \
         bioplib/FreeStringList.o  \
         bioplib/padterm.o         \
         bioplib/chindex.o         \
         bioplib/fsscanf.o         \
         bioplib/ParseRes.o        \
         bioplib/ReadPDB.o         \
         bioplib/WritePDB.o        \
         bioplib/IndexPDB.o        \
         bioplib/OpenStdFiles.o    \
         bioplib/FindResidue.o     \
         bioplib/FindResidueSpec.o \
         bioplib/FindNextResidue.o \
         bioplib/CopyPDB.o         \
         bioplib/StoreString.o     \
         bioplib/OpenFile.o        \
         bioplib/openorpipe.o      \
         bioplib/phi.o             \
         bioplib/array2.o          


cellulosetorsions : cellulosetorsions.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f cellulosetorsions.o $(LFILES)

distclean : clean
	\rm -f cellulosetorsions
