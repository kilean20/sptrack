all: main
#custom widgets
custom.o: custom.cpp custom.h
	g++ -c `wx-config --cxxflags` -I/home/hchao/usr/local/include/plplot custom.cpp

gui.o: gui.cpp gui.h custom.h
	g++ -c `wx-config --cxxflags` -I/home/hchao/usr/local/include/plplot gui.cpp

mainframe.o: mainframe.cpp mainframe.h gui.h
	g++ -c `wx-config --cxxflags` -I/home/hchao/usr/local/include/plplot mainframe.cpp

main.o: main.cpp mainframe.h gui.h
	g++ -c `wx-config --cxxflags` -I/home/hchao/usr/local/include/plplot main.cpp

main: custom.o gui.o mainframe.o main.o
	g++ `wx-config --libs` -lplplotwxwidgetsd custom.o gui.o mainframe.o main.o -o viewer

