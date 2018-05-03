Purpose of this program: 
	Convert binary outputs of the program into ASCII outputs (one file per parameter/variables), including constants
How to Compile:
	To compile the program bin2txt_params, you need to use the following command:

	g++ -I $EIGEN3_INCLUDE_DIR -I headers/  -lgsl -lgslcblas -lboost_iostreams -lboost_system -lboost_filesystem bin2txt_params.cpp sources/diagnostics.cpp sources/interpol.cpp -o bin2txt_params.out

	This assumes that the boost and gsl libraries are installed. This also assume that you have put the headers of the main program in 'headers/' and the sources in 'sources/' (subdirectories of the bin2txt_params.cpp file).

How to execute:
	./bin2txt_params.out
	The required arguments will be prompted on the screen
