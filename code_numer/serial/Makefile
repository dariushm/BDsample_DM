

CFLAGS2	     = -O3 -I /usr/include -I /Users/mjohn218/c2013.dir/rdiff/partition_code/include
LIBS2 	     = /usr/lib/libgsl.a /usr/lib/libgslcblas.a 
OBJS	     = ../code_RD/rand_gsl.o ../code_RD/md_timer.o ../code_RD/vector_rot_calls.o ../code_RD/Faddeeva.o

GC           = g++
LD           = $(CC)



survive_prob:	survive_prob.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob survive_prob.cpp $(LIBS2) $(OBJS) 
survive_prob_rate:	survive_prob_rate.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rate survive_prob_rate.cpp $(LIBS2) $(OBJS) 

survive_prob_rot:	survive_prob_rot.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rot survive_prob_rot.cpp $(LIBS2) $(OBJS) 
survive_prob_rot_thet:	survive_prob_rot_thet.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rot_thet survive_prob_rot_thet.cpp $(LIBS2) $(OBJS) 

survive_prob_rot_orient:	survive_prob_rot_orient.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rot_orient survive_prob_rot_orient.cpp $(LIBS2) $(OBJS) 

survive_prob_rotfull_orient:	survive_prob_rotfull_orient.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rotfull_orient survive_prob_rotfull_orient.cpp $(LIBS2) $(OBJS) 
survive_prob_rotfull2_orient:	survive_prob_rotfull2_orient.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rotfull2_orient survive_prob_rotfull2_orient.cpp $(LIBS2) $(OBJS) 

survive_prob_rot2_thet:	survive_prob_rot2_thet.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rot2_thet survive_prob_rot2_thet.cpp $(LIBS2) $(OBJS) 

survive_prob_rot2_orient:	survive_prob_rot2_orient.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rot2_orient survive_prob_rot2_orient.cpp $(LIBS2) $(OBJS) 

survive_prob_rot2_orient2:	survive_prob_rot2_orient2.cpp $(LIBS2) $(OBJS) 
		$(GC) $(CFLAGS2)  -o survive_prob_rot2_orient2 survive_prob_rot2_orient2.cpp $(LIBS2) $(OBJS) 
