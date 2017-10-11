#
# Makefile
#

CFLAGS		= -O3
MFLAGS		= -lm

TARGETS		= xyza2pipe ucsf2pipe nv2pipe xeasy2pipe azara2pipe vnmr2pipe xwnmr2pipe\
		pipe2xyza pipe2ucsf pipe2nv pipe2xeasy pipe2azara\
		pipe2proj add2pipe adducsf2pipe addnv2pipe addxeasy2pipe addazara2pipe addvnmr2pipe addxwnmr2pipe\
		defl2pipe

OBJECTS_MATH	= libMath.o libMatrix.o

OBJECTS_C	= initpars.o checklabel.o cnvhdr.o libString.o $(OBJECTS_MATH)

OBJECTS_XP	= checkxyza.o openxyza.o pushxyza.o

OBJECTS_PX	= checkpipe.o checkxyza.o openpipe.o pullxyza.o

OBJECTS_UP	= checkucsf.o openucsf.o pushucsf.o

OBJECTS_PU	= checkpipe.o checkucsf.o openpipe.o pullucsf.o

OBJECTS_NP	= checknv.o opennv.o pushnv.o

OBJECTS_PN	= checkpipe.o checknv.o openpipe.o pullnv.o

OBJECTS_EP	= checkxeasy.o openxeasy.o pushxeasy.o xeasy2float.o

OBJECTS_PE	= checkpipe.o checkxeasy.o openpipe.o pullxeasy.o xeasy2float.o

OBJECTS_AP	= checkazara.o openazara.o pushazara.o

OBJECTS_PA	= checkpipe.o checkazara.o openpipe.o pullazara.o

OBJECTS_VP	= checkvnmr.o vendorpar.o openvnmr.o pushvnmr.o

OBJECTS_BP	= checkxwnmr.o vendorpar.o openxwnmr.o pushxwnmr.o

OBJECTS_PJ	= checkpipe.o checkxyza.o openpipe.o pullproj.o

OBJECTS_DXP	= checkxyza.o openxyza.o pushadd.o

OBJECTS_DUP	= checkucsf.o openucsf.o pushadducsf.o

OBJECTS_DNP	= checknv.o opennv.o pushaddnv.o

OBJECTS_DEP	= checkxeasy.o openxeasy.o pushaddxeasy.o xeasy2float.o

OBJECTS_DAP	= checkazara.o openazara.o pushaddazara.o

OBJECTS_DVP	= checkvnmr.o vendorpar.o openvnmr.o pushaddvnmr.o

OBJECTS_DBP	= checkxwnmr.o vendorpar.o openxwnmr.o pushaddxwnmr.o

OBJECTS_DFL	= checkdefl.o openxyza.o pushxyza.o

BIN		= ./bin/

ADD2PIPE_ALIAS	= addxyza2pipe

all: $(TARGETS)
	mkdir -p $(BIN)
	cp -f $(TARGETS) $(BIN)
	@if [ ! -f $(BIN)$(ADD2PIPE_ALIAS) ]; then (cd $(BIN); ln -s add2pipe $(ADD2PIPE_ALIAS)); fi

clean:
	rm -f *.o
	rm -f $(TARGETS) $(ADD2PIPE_ALIAS)
	rm -f $(addprefix $(BIN), $(TARGETS) addxyza2pipe)

$(OBJECTS_MATH):
	$(CC) $*.c -c -o $@ $(CFLAGS)

.o:
	$(CC) $< -c -o $@ $(CFLAGS)

xyza2pipe: $(OBJECTS_C) $(OBJECTS_XP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

pipe2xyza: $(OBJECTS_C) $(OBJECTS_PX)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

ucsf2pipe: $(OBJECTS_C) $(OBJECTS_UP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

pipe2ucsf: $(OBJECTS_C) $(OBJECTS_PU)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

nv2pipe: $(OBJECTS_C) $(OBJECTS_NP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

pipe2nv: $(OBJECTS_C) $(OBJECTS_PN)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

xeasy2pipe: $(OBJECTS_C) $(OBJECTS_EP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

pipe2xeasy: $(OBJECTS_C) $(OBJECTS_PE)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

azara2pipe: $(OBJECTS_C) $(OBJECTS_AP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

pipe2azara: $(OBJECTS_C) $(OBJECTS_PA)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

vnmr2pipe: $(OBJECTS_C) $(OBJECTS_VP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

xwnmr2pipe: $(OBJECTS_C) $(OBJECTS_BP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

pipe2proj: $(OBJECTS_C) $(OBJECTS_PJ)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

add2pipe: $(OBJECTS_C) $(OBJECTS_DXP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)
	test ! -f $(ADD2PIPE_ALIAS) && ln -s $@ $(ADD2PIPE_ALIAS)

adducsf2pipe: $(OBJECTS_C) $(OBJECTS_DUP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

addnv2pipe: $(OBJECTS_C) $(OBJECTS_DNP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

addxeasy2pipe: $(OBJECTS_C) $(OBJECTS_DEP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

addazara2pipe: $(OBJECTS_C) $(OBJECTS_DAP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

addvnmr2pipe: $(OBJECTS_C) $(OBJECTS_DVP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

addxwnmr2pipe: $(OBJECTS_C) $(OBJECTS_DBP)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)

defl2pipe: $(OBJECTS_C) $(OBJECTS_DFL)
	$(CC) $@.c $^ -o $@ \
	$(CFLAGS) $(MFLAGS)
