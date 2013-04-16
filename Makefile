# change to your local directories!
PD_APP_DIR = /Applications/Pd-extended.app/Contents/Resources
GEM_DIR = /Users/hans/auto-build/pd-extended/externals/Gem

CPPFLAGS = -I$(GEM_DIR)/src -I$(PD_APP_DIR)/include

#linux doesnt work yet
UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
 CPPFLAGS += 
 CFLAGS = -g -O2 -fPIC -freg-struct-return -Os -falign-loops=32 -falign-functions=32 -falign-jumps=32 -funroll-loops -ffast-math -mmmx -fpascal-strings 
 LDFLAGS = -shared -rdynamic
 LIBS = -lm
 EXTENSION = pd_linux
 USER_EXTERNALS=$(HOME)/pd-externals
endif
ifeq ($(UNAME),Darwin)
 CPPFLAGS += 
 CFLAGS = -g -fast -msse3 -arch i386
 LDFLAGS = -bundle -bundle_loader $(PD_APP_DIR)/bin/pd -undefined dynamic_lookup -arch i386 -mmacosx-version-min=10.6
 LIBS = -lm
 EXTENSION = pd_darwin
 USER_EXTERNALS=$(HOME)/Library/Pd
endif

.SUFFIXES = $(EXTENSION)

SOURCES = dbap3d
# dbap2d

all:
	gcc $(CPPFLAGS) $(CFLAGS) -o $(SOURCES).o -c $(SOURCES).c
	gcc -o $(SOURCES).$(EXTENSION) $(LDFLAGS) $(SOURCES).o $(LIBS)
	rm -fr ./*.o

deploy:
	install -d $(USER_EXTERNALS)/$(SOURCES)
	install -p *.pd $(USER_EXTERNALS)/$(SOURCES)
	install -p $(SOURCES).$(EXTENSION) $(USER_EXTERNALS)/$(SOURCES)

clean:
	rm -f $(SOURCES)*.o
	rm -f $(SOURCES)*.$(EXTENSION)

distro: clean all
	rm *.o
