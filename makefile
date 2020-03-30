CGALCONFIG=
CGALINCLUDE=/usr/include/CGAL
CGALLIB=/usr/lib/x86_64-linux-gnu
BOOSTINCLUDE=/usr/include/boost
BOOSTLIB=/usr/lib/x86_64-linux-gnu

all: smatch mtgen compgen

smatch: smatch.cpp SES.o contextBall.o util.o gridSES.o memusage.o MyDockinfo.o
	g++ -I$(BOOSTINCLUDE) -I$(CGALCONFIG) -I$(CGALINCLUDE) -O3 -Wall -c -fmessage-length=0 smatch.cpp
	g++ -static -L$(BOOSTLIB) -L$(CGALLIB) -I$(BOOSTINCLUDE) -I$(CGALCONFIG) -I$(CGALINCLUDE) -O3 -o smatch smatch.cpp util.o MyDockinfo.o SES.o contextBall.o gridSES.o memusage.o -lCGAL -lm -lboost_serialization -lpthread

SES.o: SES.cpp SES.h global.h
	g++ -I$(BOOSTINCLUDE) -I$(CGALCONFIG) -I$(CGALINCLUDE) -O3 -Wall -c -fmessage-length=0 SES.cpp
	
gridSES.o: gridSES.h
	g++ -I$(BOOSTINCLUDE) -I$(CGALCONFIG) -I$(CGALINCLUDE) -O3 -Wall -c -fmessage-length=0 gridSES.cpp

contextBall.o: contextBall.cpp contextBall.h util.h global.h
	g++ -I$(BOOSTINCLUDE) -I$(CGALCONFIG) -I$(CGALINCLUDE) -O3 -Wall -c -fmessage-length=0 contextBall.cpp

MyDockinfo.o: MyDockinfo.h MyDockinfo.cpp util.h
	g++ -I$(BOOSTINCLUDE) -I$(CGALCONFIG) -I$(CGALINCLUDE) -O3 -Wall -c -fmessage-length=0 MyDockinfo.cpp

util.o: util.h util.cpp
	g++ -O3 -c util.cpp

memusage.o: memusage.h memusage.cpp
	g++ -O3 -c memusage.cpp

mtgen: mtgen.cpp
	g++ -I$(BOOSTINCLUDE) -I$(CGALCONFIG) -I$(CGALINCLUDE) -O3 -Wall -c -fmessage-length=0 mtgen.cpp
	g++ -static -L$(BOOSTLIB) -L$(CGALLIB) -O3 -o mtgen mtgen.o -lCGAL -lm -lboost_serialization -lpthread

clean:
	rm smatch mtgen compgen *.o *.gch
