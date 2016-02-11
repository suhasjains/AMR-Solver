NAME 	:= myoctree
OBDIR	:= ../OBJ/
RUNDIR	:= ../RUN/

OBJECTS =   $(OBDIR)myoctree.o $(OBDIR)octree.o $(OBDIR)block.o $(OBDIR)output.o 

HEADER 	:= octree.h block.h output.h 

.DEFAULT:

INCPATH	=	-I ~/ 

LIBPATH	= 	-L ~/

LIBS	=	-lm  

CXX	:=	mpic++
CXXFLAGS:= 	-g  $(INCPATH) -Wall -O3
#FLAGS	:=	-g -pg $(INCPATH) -Wall -fopenmp -O3
#LFLAGS	:=	$(LIBPATH) $(LIBS) -lm -lgsl -lblas
LFLAGS	:=	$(LIBPATH) $(LIBS) 

####################################################################

$(NAME): $(OBJECTS) $(RUNDIR)
	$(CXX) -o $(NAME) $(OBJECTS) $(FLAGS) $(LFLAGS)
	@mv $(NAME) $(RUNDIR) 
#	@cp iwrite $(RUNDIR)

$(RUNDIR): 
	@test -d $(RUNDIR) || mkdir $(RUNDIR)

$(OBDIR)myoctree.o: myoctree.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)myoctree.o myoctree.cpp

$(OBDIR)octree.o: octree.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)octree.o octree.cpp

$(OBDIR)block.o: block.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)block.o block.cpp

$(OBDIR)output.o: output.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)output.o output.cpp
copy:
	cp -f clean $(RUNDIR)/.	
clean: 		
	$(RM) *~
	$(RM) $(OBDIR)*.o
