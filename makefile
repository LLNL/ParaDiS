#-------------------------------------------------------------------------------------------
#
#    makefile:  This is the primary make file controlling the build
#               of the ParaDiS parallel code and associated utilities.
#
#    Usage:
#        make           build paradis executable and support utilities
#        make clean     remove executable and object files
#        make depend    update makefile dependencies
#
#-------------------------------------------------------------------------------------------

BINDIR = ./bin
OBJDIR = ./obj

#-------------------------------------------------------------------------------------------

.PHONY: all
all: $(OBJDIR) $(BINDIR)
	( cd src   ; $(MAKE) -j; )
	( cd utils ; $(MAKE) -j; )

.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(BINDIR) .DS_Store */.DS_Store
	( cd src   ; $(MAKE) clean ; )
	( cd utils ; $(MAKE) clean ; )
	( cd tests ; $(MAKE) clean ; )
	( cd python; $(MAKE) clean ; )

.PHONY: ext
ext:
	( cd ./ext; $(MAKE); )

.PHONY: depend
depend:
	( cd src ; $(MAKE) $@ ; cd .. ; )

.PHONY: purify
purify: $(BINDIR)
	( cd src ; $(MAKE) $@ ; cd .. ; )

.PHONY: prof
prof: $(BINDIR)
	( cd src ; $(MAKE) $@ ; cd .. ; )

.PHONY: rebuild
rebuild: clean $(OBJDIR) $(BINDIR)
	( cd src   ; $(MAKE) -j; )
	( cd utils ; $(MAKE) -j; )

.PHONY: rebuild_all
rebuild_all: clean $(OBJDIR) $(BINDIR)
	( cd ext   ; $(MAKE) clean ; $(MAKE) ; )
	( cd src   ; $(MAKE) -j; )
	( cd utils ; $(MAKE) -j; )

$(OBJDIR) :
	mkdir -p $(OBJDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)
