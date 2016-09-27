
CFLAGS	     = -O3 -I$(MKLROOT)/include
CC           = icpc
LD           = $(CC)

CFLAGS2	     = -O3 -I /data/johnsonme/rdiff/include -I /home/johnsonme/hhpc/partition_code/include -I /data/johnsonme/gsl/include
##LIBS2 	     = /usr/lib64/libgsl.so /usr/lib64/libgslcblas.so -lrt  
LIBS2 	     = /data/johnsonme/gsl/lib/libgsl.a
OBJS	     = rand_gsl.o /data/johnsonme/gsl/lib/libgsl.a /data/johnsonme/clathrin/code/md_timer.o /home/johnsonme/hhpc/partition_code/code_RD/vector_rot_calls.o /home/johnsonme/hhpc/partition_code/code_RD/Faddeeva.o
GC	     = g++



pde_gf_radial_mkl:	pde_gf_radial_mkl.o
		$(LD) -o $@ pde_gf_radial_mkl.o $(LIBS)


survive_prob:	survive_prob.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob survive_prob.cpp $(LIBS2) $(OBJS) 

survive_prob_rot:	survive_prob_rot.cpp $(OBJS) $(LIBS2)
		$(GC) $(CFLAGS2)  -o survive_prob_rot survive_prob_rot.cpp $(OBJS) $(LIBS2)
survive_prob_rate:	survive_prob_rate.cpp $(OBJS) $(LIBS2)
		$(GC) $(CFLAGS2)  -o survive_prob_rate survive_prob_rate.cpp $(OBJS) $(LIBS2)
survive_prob_rot_thet:	survive_prob_rot_thet.cpp $(OBJS) $(LIBS2)
		$(GC) $(CFLAGS2)  -o survive_prob_rot_thet survive_prob_rot_thet.cpp $(OBJS) $(LIBS2)

survive_prob_rot_orient:	survive_prob_rot_orient.cpp $(OBJS) $(LIBS2)
		$(GC) $(CFLAGS2)  -o survive_prob_rot_orient survive_prob_rot_orient.cpp $(OBJS) $(LIBS2) 

survive_prob_rot2_thet:	survive_prob_rot2_thet.cpp $(OBJS) $(LIBS2)
		$(GC) $(CFLAGS2)  -o survive_prob_rot2_thet survive_prob_rot2_thet.cpp $(OBJS) $(LIBS2)

survive_prob_rot2_orient:	survive_prob_rot2_orient.cpp $(OBJS) $(LIBS2)
		$(GC) $(CFLAGS2)  -o survive_prob_rot2_orient survive_prob_rot2_orient.cpp $(OBJS) $(LIBS2)
survive_prob_rotfull_orient:	survive_prob_rotfull_orient.cpp $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rotfull_orient survive_prob_rotfull_orient.cpp $(OBJS) 

survive_prob_rot2_orient2:	survive_prob_rot2_orient2.cpp $(OBJS) $(LIBS2)
		$(GC) $(CFLAGS2)  -o survive_prob_rot2_orient2 survive_prob_rot2_orient2.cpp $(OBJS) $(LIBS2)
