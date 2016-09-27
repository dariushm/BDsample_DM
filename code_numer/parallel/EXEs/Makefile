HOME	= /Users/mjohn218
FPRDIR= /home/likewise-open/WIN/oyogurt1/workspace/PBC/2dSIMFINAL/fpr_code-submaster1/FPR_CELL
##$(HOME)/github_dirs/fpr_code/FPR_CELL

LD           = $(CC)
MPCC         = mpicxx
CFLAGS	= -O3 -I ../include -I $(FPRDIR)/include
CFLAGS += $(shell gsl-config --cflags)
LIBS	= $(shell gsl-config --libs)

SDIR	= ../subroutines
OBJS	     = $(SDIR)/rand_gsl.o $(SDIR)/md_timer.o $(SDIR)/crossproduct.o $(SDIR)/rotationEuler.o $(SDIR)/rotate.o
OBJS1	     = $(SDIR)/angle_calc.o $(SDIR)/GaussV.o $(SDIR)/calc_Rmatrix.o

CC	     = g++



par_survive_rot2:	par_survive_rot2.cpp  $(OBJS) 
		$(MPCC) $(CFLAGS)  -o par_survive_rot2 par_survive_rot2.cpp $(OBJS) $(LIBS)
par_rate_rotleg:	par_rate_rotleg.cpp $(OBJS)
		$(MPCC) $(CFLAGS)  -o par_rate_rotleg par_rate_rotleg.cpp $(OBJS) $(LIBS)

par_survive_vr:	par_survive_vr.cpp $(OBJS)
		$(MPCC) $(CFLAGS)  -o par_survive_vr par_survive_vr.cpp $(OBJS) $(LIBS)
par_survive_vr_cut:	par_survive_vr_cut.cpp $(OBJS)
		$(MPCC) $(CFLAGS)  -o par_survive_vr_cut par_survive_vr_cut.cpp $(OBJS) $(LIBS)
par_survive_vrdt:	par_survive_vrdt.cpp $(OBJS)
		$(MPCC) $(CFLAGS)  -o par_survive_vrdt par_survive_vrdt.cpp $(OBJS) $(LIBS)

##below is just temporary, set to do costheta=0,0.
par_survive_rot2_0:	par_survive_rot2_0.cpp 
		$(MPCC) $(CFLAGS)  -o par_survive_rot2_0 par_survive_rot2_0.cpp $(OBJS) $(LIBS)

par_survive:	par_survive.cpp 
		$(MPCC) $(CFLAGS)  -o par_survive par_survive.cpp $(OBJS) $(LIBS)
par_survive_mem:	par_survive_mem.cpp 
		$(MPCC) $(CFLAGS)  -o par_survive_mem par_survive_mem.cpp $(OBJS) $(LIBS)

par_survive_rotfull_orient:	par_survive_rotfull_orient.cpp 
		$(MPCC) $(CFLAGS)  -o par_survive_rotfull_orient par_survive_rotfull_orient.cpp $(OBJS) $(LIBS)

par_survive_rotcom_orient:	par_survive_rotcom_orient.cpp 
		$(MPCC) $(CFLAGS)  -o par_survive_rotcom_orient par_survive_rotcom_orient.cpp $(OBJS1) $(OBJS2) $(LIBS)
par_survive2_rotcom_orient:	par_survive2_rotcom_orient.cpp $(OBJS1)
		$(MPCC) $(CFLAGS)  -o par_survive2_rotcom_orient par_survive2_rotcom_orient.cpp $(OBJS1) $(OBJS2) $(LIBS)
par_survive3_rotcom_orient:	par_survive3_rotcom_orient.cpp $(OBJS1)
		$(MPCC) $(CFLAGS)  -o par_survive3_rotcom_orient par_survive3_rotcom_orient.cpp $(OBJS1) $(OBJS2) $(LIBS)
par_sadapt_rotcom_orient:	par_sadapt_rotcom_orient.cpp $(OBJS1)
		$(MPCC) $(CFLAGS)  -o par_sadapt_rotcom_orient par_sadapt_rotcom_orient.cpp $(OBJS1) $(OBJS2) $(LIBS)
par_survive4_rotcom_orient:	par_survive4_rotcom_orient.cpp $(OBJS1)
		$(MPCC) $(CFLAGS)  -o par_survive4_rotcom_orient par_survive4_rotcom_orient.cpp $(OBJS1) $(OBJS2) $(LIBS)

par_rate_rotcom_orient:	par_rate_rotcom_orient.cpp $(OBJS1) $(OBJS)
		$(MPCC) $(CFLAGS)  -o par_rate_rotcom_orient par_rate_rotcom_orient.cpp $(OBJS) $(OBJS1) $(LIBS)
par_rate2_rotcom_orient:	par_rate2_rotcom_orient.cpp $(OBJS1)
		$(MPCC) $(CFLAGS)  -o par_rate2_rotcom_orient par_rate2_rotcom_orient.cpp $(OBJS1) $(OBJS2) $(LIBS)
par_rate3_rotcom_orient:	par_rate3_rotcom_orient.cpp $(OBJS1)
		$(MPCC) $(CFLAGS)  -o par_rate3_rotcom_orient par_rate3_rotcom_orient.cpp $(OBJS1) $(OBJS2) $(LIBS)

par_survive_rotfull_orient_rename:	par_survive_rotfull_orient_rename.cpp 
		$(MPCC) $(CFLAGS)  -o par_survive_rotfull_orient_rename par_survive_rotfull_orient_rename.cpp $(OBJS) $(OBJS2) $(LIBS) 

par_survive_rotfull_orient_phi12:	par_survive_rotfull_orient_phi12.cpp 
		$(MPCC) $(CFLAGS)  -o par_survive_rotfull_orient_phi12 par_survive_rotfull_orient_phi12.cpp $(OBJS) $(OBJS2) $(LIBS) 


$(SDIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $<"
	$(CC) $(CFLAGS) -c $< -o $@  
	@echo "------------"
