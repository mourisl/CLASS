BAM_INCLUDE = -I./samtools-0.1.19/
LP_SOLVE_INCLUDE = -I./lp_solve_5.5
BAM_LIB = -L./samtools-0.1.19/ -L/cm/shared/apps/zlib/1.2.8/lib/
all: Reads.o SplicesInformation.o SolveRegion.o FindRegions.o TranscriptDecider.o BitTable.o Evidence.o TTS.o main.cpp 
	g++ -O -o class main.cpp ${BAM_LIB} ${LP_SOLVE_INCLUDE} ${BAM_INCLUDE} Reads.o SplicesInformation.o \
		SolveRegion.o FindRegions.o TranscriptDecider.o Evidence.o BitTable.o TTS.o liblpsolve55.a -lbam -lz -lm -lpthread -ldl
	g++ -O -o junc FindJunction.cpp  ${BAM_LIB} ${BAM_INCLUDE} -lbam -lz -lm -lpthread
	g++ -O -o clnb CleanBam.cpp  ${BAM_LIB} ${BAM_INCLUDE} -lbam -lz -lm -lpthread
	g++ -O -o grader grader.cpp
	g++ -O -o addXS AddXS.cpp
SolveRegion.o: SolveRegion.cpp SolveRegion.h 
	g++ ${BAM_INCLUDE} -O -c SolveRegion.cpp ${LP_SOLVE_INCLUDE} 
Reads.o: Reads.cpp Reads.h
	g++ -O -c Reads.cpp ${BAM_INCLUDE}  
FindRegions.o: FindRegions.cpp FindRegions.h
	g++ -O -c FindRegions.cpp ${LP_SOLVE_INCLUDE} ${BAM_INCLUDE} 
SplicesInformation.o: SplicesInformation.cpp SplicesInformation.h
	g++ -O -c SplicesInformation.cpp 
Evidence.o: Evidence.cpp Evidence.h
	g++ -O -c Evidence.cpp ${LP_SOLVE_INCLUDE} ${BAM_INCLUDE}  
TranscriptDecider.o: TranscriptDecider.cpp TranscriptDecider.h
	g++ -O -c TranscriptDecider.cpp ${LP_SOLVE_INCLUDE} ${BAM_INCLUDE}  
TTS.o: TTS.cpp TTS.h
	g++ -O -c TTS.cpp 
BitTable.o: BitTable.cpp BitTable.h
	g++ -O -c BitTable.cpp BitTable.h
clean:
	rm *.o *.gch ./class ./junc ./clnb ./grader ./addXS
