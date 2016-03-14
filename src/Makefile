NAME 	:= amrsolver
OBDIR	:= ../obj/
RUNDIR	:= ../output/

OBJECTS =   $(OBDIR)amrsolver.o $(OBDIR)octreegrid.o $(OBDIR)adapt.o $(OBDIR)octree.o $(OBDIR)block.o $(OBDIR)vecfield.o $(OBDIR)field.o $(OBDIR)boundary.o $(OBDIR)ghost.o $(OBDIR)vtk.o $(OBDIR)output.o $(OBDIR)input.o 

HEADER 	:= octreegrid.h adapt.h octree.h block.h field.h vecfield.h vtk.h output.h boundary.h ghost.h direction.h 

.DEFAULT:

INCPATH	=	-I ~/ 

LIBPATH	= 	-L ~/

LIBS	=	-lm  

CXX	:=	mpic++
CXXFLAGS:= 	-g  -Wall -O3 -fno-elide-constructors
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

$(OBDIR)amrsolver.o: amrsolver.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)amrsolver.o amrsolver.cpp

$(OBDIR)octreegrid.o: octreegrid.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)octreegrid.o octreegrid.cpp

$(OBDIR)adapt.o: adapt.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)adapt.o adapt.cpp

$(OBDIR)octree.o: octree.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)octree.o octree.cpp

$(OBDIR)field.o: field.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)field.o field.cpp

$(OBDIR)vecfield.o: vecfield.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)vecfield.o vecfield.cpp

$(OBDIR)block.o: block.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)block.o block.cpp

$(OBDIR)boundary.o: boundary.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)boundary.o boundary.cpp

$(OBDIR)ghost.o: ghost.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)ghost.o ghost.cpp

$(OBDIR)input.o: input.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)input.o input.cpp

$(OBDIR)vtk.o: vtk.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)vtk.o vtk.cpp

$(OBDIR)output.o: output.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)output.o output.cpp
copy:
	cp -f clean $(RUNDIR)/.	
clean: 		
	$(RM) *~
	$(RM) $(OBDIR)*.o
