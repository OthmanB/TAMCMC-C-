### Purpose of getmodel ###
	Read a simple input file that contains length and a list of vector of parameters. From this, it builds a model, provided that the user gave:
			- the model name
			- the data filename, in the same format as when given to TAMCMC
	The model is written in an output ASCII file. That output can be user-specified, or a default one


### Purpose of bin2txt ###
	Convert binary outputs of the program into ASCII outputs (one file per parameter/variables), including constants
