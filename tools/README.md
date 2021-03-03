### Purpose of getmodel ###
	Read a simple input file that contains length and a list of vector of parameters. From this, it builds a model, provided that the user gave:
			- the model name
			- the data filename, in the same format as when given to TAMCMC
	The model is written in an output ASCII file. That output can be user-specified, or a default one


### Purpose of bin2txt ###
	Convert binary outputs of the program into ASCII outputs (one file per parameter/variables), including constants

### Purpose of read_stats ###
    Convert binary files with the Likelihood/Prior/Posterior into ASCII outputs.
    
    
### Purpose of a2_nl.py ###
    To be used only for models that handle a2 coefficient using a polynomial function. This program translates the a2 polynomial coeficients extracted using data_extract_IDL.zip into a2(nu_nl) with a proper handling of the error propagation. It makes a jpg plot of the this.

