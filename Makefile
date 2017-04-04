all:
	g++ neur27.cxx rk.cxx -o exe
	./exe input27

plot:
	matlab -nodisplay -nosplash -r "plotit;quit"
