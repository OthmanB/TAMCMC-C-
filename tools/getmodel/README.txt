Purpose of this program: 
	Read a simple input file that contains length and a list of vector of parameters. From this, it builds a model, provided that the user gave:
			- the model name
			- the data filename, in the same format as when given to TAMCMC
	The model is written in an output ASCII file. That output can be user-specified, or a default one

How to Compile:
	To compile the program getmodel, you need to use the following command:

	g++ -I $EIGEN3_INCLUDE_DIR -I headers/ -lgsl -lgslcblas -lboost_iostreams -lboost_system -lboost_filesystem sources/models.cpp sources/string_handler.cpp sources/config.cpp sources/noise_models.cpp sources/function_rot.cpp sources/build_lorentzian.cpp sources/interpol.cpp sources/matrices.cpp sources/io_ms_global.cpp getmodel.cpp -o getmodel.out


	This assumes that the boost and gsl libraries are installed and that you have proper links to the sources and header directories of the main program.

	NOTE: The compilation will fail if you try to compile from another directory than the default directory. This is also true for the execution!
How to execute:
	./gemodel.out
	The required arguments will be prompted on the screen
