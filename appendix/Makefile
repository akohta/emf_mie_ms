CC      =gcc
CFLAGS  =-O3 -Wall
LDFLAGS =-lm -lgsl -lgslcblas
MPFLAGS =-fopenmp
SRCDIR1 =../exmie_src
SRCDIR2 =../mfb_src
SRCDIR3 =../com_src
OBJDIR  =../obj
TARGET2 =radiation_force.out
TRGSRC2 =radiation_force.c

SRCS1=$(wildcard $(SRCDIR1)/*.c)
OBJS1=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS1))))
HEAD1=$(wildcard $(SRCDIR1)/*.h)

SRCS2=$(wildcard $(SRCDIR2)/*.c)
OBJS2=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS2))))
HEAD2=$(wildcard $(SRCDIR2)/*.h)

SRCS3=$(wildcard $(SRCDIR3)/*.c)
OBJS3=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS3))))
HEAD3=$(wildcard $(SRCDIR3)/*.h)


TRGOBJ2=$(filter-out $(OBJDIR)/mie_ms_solver.o,$(OBJS1)) $(OBJS2) $(OBJS3) $(patsubst %.c,%.o,$(TRGSRC2))

all : directories $(TARGET2)

directories:
	@mkdir -p $(OBJDIR)

$(TARGET2) : $(TRGOBJ2)
	$(CC) $(LDFLAGS) $(MPFLAGS) -o $@ $^

$(OBJDIR)/%.o : $(SRCDIR1)/%.c
	$(CC) $(CFLAGS) $(MPFLAGS) -I$(SRCDIR2) -I$(SRCDIR3) -c $< -o $@
	
$(OBJDIR)/%.o : $(SRCDIR2)/%.c
	$(CC) $(CFLAGS) -I$(SRCDIR3)-c $< -o $@

$(OBJDIR)/%.o : $(SRCDIR3)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

.c.o :
	$(CC) $(CFLAGS) $(MPFLAGS) -I$(SRCDIR1) -I$(SRCDIR2) -I$(SRCDIR3) -c $<
	
clean:
	@rm -rf $(TARGET2) ./*.o
	
$(OBJS1) : $(HEAD1) $(HEAD2) $(HEAD3)
$(OBJS2) : $(HEAD2) $(HEAD3)
$(OBJD3) : $(HEAD3)
