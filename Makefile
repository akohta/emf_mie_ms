# Makefile for intel C compiler  
CC      =gcc
CFLAGS  =-O3 -Wall
LDFLAGS =-lm -lgsl -lgslcblas
MPFLAGS =-fopenmp
SRCDIR1 =exmie_src
SRCDIR2 =mfb_src
OBJDIR  =obj
TARGET1 =mie_ms_solver
TRGSRC1 =mie_ms_solver.c
TARGET2 =example1.out
TRGSRC2 =example1.c
TARGET3 =example2.out
TRGSRC3 =example2.c

SRCS1=$(wildcard $(SRCDIR1)/*.c)
OBJS1=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS1)) ))
HEAD1=$(wildcard $(SRCDIR1)/*.h)

SRCS2=$(wildcard $(SRCDIR2)/*.c)
OBJS2=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS2)) ))
HEAD2=$(wildcard $(SRCDIR2)/*.h)

TRGOBJ1=$(OBJS1) $(OBJS2) 
TRGOBJ2=$(filter-out $(OBJDIR)/$(TARGET1).o,$(OBJS1)) $(OBJS2) $(patsubst %.c,%.o,$(TRGSRC2))
TRGOBJ3=$(filter-out $(OBJDIR)/$(TARGET1).o,$(OBJS1)) $(OBJS2) $(patsubst %.c,%.o,$(TRGSRC3))

all : directories $(TARGET1) $(TARGET2) $(TARGET3)

directories:
	@mkdir -p $(OBJDIR)

$(TARGET1) : $(TRGOBJ1)  
	$(CC) $(LDFLAGS) $(MPFLAGS) -o $@ $^

$(TARGET2) : $(TRGOBJ2)  
	$(CC) $(LDFLAGS) $(MPFLAGS) -o $@ $^

$(TARGET3) : $(TRGOBJ3)  
	$(CC) $(LDFLAGS) $(MPFLAGS) -o $@ $^

$(OBJDIR)/%.o : $(SRCDIR1)/%.c
	$(CC) $(CFLAGS) $(MPFLAGS) -I $(SRCDIR2) -c $< -o $@
	
$(OBJDIR)/%.o : $(SRCDIR2)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
	
.c.o :
	$(CC) $(CFLAGS) $(MPFLAGS) -I$(SRCDIR1) -I$(SRCDIR2) -c $<
	
clean:
	@rm -rf $(TARGET1) $(TARGET2) $(TARGET3) $(OBJDIR) ./*.o
	
$(OBJS1) : $(HEAD1) $(HEAD2)
$(OBJS2) : $(HEAD2)
