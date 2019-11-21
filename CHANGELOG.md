# Version history #

### v1.3.3-dev New model and improvements ###
	  * Add two models with fit of the Width from Appourchaux 2012 instead of individual widths. 
	    This is made in a transparent way, using the input widths for defining the guesses. Thus, no need to change drastically the structure of .model files
	  	Priors on the parameter of the model are hard-coded at the moment. See io_ms_global.cpp for their value (typically 10-20% of the expected values) [DONE] [NEED TEST ON DALMA]
	  * When starting a new process, a 'version.txt' file is writtent in the object directory. This file gives the code version used for the processing [DONE] [TESTED]
	  * Adding a hard-coded limit on a3/a1 < 0.2. This to avoid very unphysical solutions (that lead to a star rotating in the opposite direction in the pole [DONE] [TESTED]
	  * Redesigning the way likelihoods, models, priors and prime priors are handled: Instead of having them hard coded, they are now defined along with their switch case into *.list 
	    files within Config/default. This significantly simplifies the implementation of new models.
	    
### v1.3.2-dev Improvements ###
      * Adding the possibility to fit amplitudes instead of Height by specifying  fit_squareAmplitude_instead_Height   [bool value]   into the .model file [DONE] [TESTED]
      * Add the possibility to change the priors for Height/Amplitudes/Width in the .model file [DONE] [NEED THOROUGH COMPARATIVE TESTING WITH EARLIER STABLE VERSION] 
        Note that we are permissive regarding the denomination Height or Amplitude. However, an explicit fit of squaredAmplitude could be requested. See example of .model file
      * Add the possibility to change frequency prior between GUG or Uniform. [DONE] [NEED THOROUGH COMPARATIVE TESTING WITH EARLIER STABLE VERSION] 
        For GUG the syntax is
      			Frequency     GUG     -1    -1    -1     [sigma1]    [sigma2]       (-1 are replaced by values of the table of frequencies)
      	For Uniform the syntax is
      			Frequency     Uniform  -1  -1   -1      (-1 are replaced by values of the table of frequencies)
      * Reorganising io_ms_global using new functions: initialise_params(), fill_param() and add_params() 
        This reduced the size of this program significantly by removing redoundant commands
        This allows much easier generation of a model. Note that the function could be used to create any kind of vector of parameters
      * Minor esthetic improvements in the text outputs
      * Few minor typos and fixes into the .md files [DONE]
      * Starting a basic documentation in tex [DEV] [POSTPONED: Replace by Wifi on Github]
        
### v1.3.1 Minor improvements/Bug fix ###
	* Correcting a typo in getmodel.cpp that prevented the compilation of getmodel tool
	* Further cleaning of useless lines, in particular in outputs.cpp
	* Removal of a debug text shown when reading the header of the parameters
	* Adding model_MS_Global_a1etaa3_Classic which allows to directly fit a1 and inclination (rather than sqrt(a1).cosi and sqrt(a1).sini)
	* Adding adaptive scheme for MS_Global_a1eta3: If the user provide splitting_a1 and inclination instead of the expected parameters (sqrt(a1)cosi, sqrt(a1)sini), then 
	  the program convert and adjust the inputs so that they are compatible with the model
	* Removal of some lines that were added in v1.3 and that stop the program if using 4-lines format for the s2 inputs (instead of 3 in simulations)
	* Adding a failsafe if the user request to read a y-data column that does not exists
      * FUTURE:
  		* Write comprehensive md files /doc file with compilation configuration and execution instructions 
  		* Verifying the Gaussian fit case
  		* Handling of complex space fitting
  		* Continue clean handling of errors perhaps using 'Either'? see https://hackernoon.com/error-handling-in-c-or-why-you-should-use-eithers-in-favor-of-exceptions-and-error-codes-f0640912eb45
  				* Create a dedicated message.cpp / message.h that handles messages and errors (move Diagnostics::file_error in there)		 	
  				* Better structure for output.cpp ... use header to save constants


###  v.1.3 Major Improvements ###
	* Detection of openmp and option to deactivate it in order to compile with clang that is not compatible with openmp
	* Detection of GSL and option to deactivate it if this library is not available on the considered machine
	* Detection of GNUPLOT. It is mandatory
	* Option to compile on DALMA (unified cmakefile)
	* Option to show the version of the program/compiler added.
	* Variables start_index_ids and last_index_ids are now obselete. To select process to be run, these are now passed as argument of the executable
	* Better handling of io_MS_Global 
		* In the config_default.cfg, only prior_fct_name is now required to specify which model should be used. 
		When prior_fct_name=io_ms_Global (capital letters matter!) the code calls io_ms_global.cpp, that read the .model file. The .model file
		must now contain a line that specify the model name (common variable model_fullname).  
		Basically, prior_fct_name can act as a switch between models in the future
		* Use of an error function for repetitive error messages 
		* The variable average_a1nl is now useless in .model files. This due to the proper handling of the model names in io_ms_global.cpp
		* The c constant that controls the Lorentzian truncation is now replaced by a variable called trunc_c that could be added in the .model file
		  If set to -1, trunc_c=10000., a very large value that is equivalent to no truncation at all.
	* Cleaner Diagnostic.cpp / config.cpp / config_preset.cpp / MALA.cpp
		* Use of an error function for file errors at opening time. 
		* Automatic handling of missing gsl (a single diagnostic.cpp) 
		* Function strsplit, strtrim, dbl_to_str, VectXd_to_Vec, str_to_dblarr, str_to_Xdarr, int_to_str are inside string_handler.cpp and not duplicated 
	* bin2txt and getmodel:
		. Complete integration within the main project, including compilation.
		. Version option support for bin2txt and getmodel 
	* getstats integration within the main project 
		. Adding function that check consistency with MS_Global models that predate 1.3.0
		. Adding function that adapt input vector if plength=10 and the we deal with an MS_Global model

### v1.2.3 Improvements ###
	* Added functionalities:
		The keyword 'average_a1nl' was added. It allows you to decide wether to fit a1(n), a1(n,l) or a1(l).
		NOTE THAT YOU DO NOT HAVE A CONTROL ON THE FIXED/FREE VARIABLES THERE. e.g. if you want to fix a1(l=1) = a11 and let free a1(l=2), you must do it by code edition
		
		The syntax should be as follow: average_a1nl     bool    [0/1]   [0/1]
        
        Ex 1: average_a1nl     bool    1    1" << std::endl;
              This uses model_MS_Global_a1etaa3_* which assumes a1(n,l) = a1
        
        Ex 2: average_a1nl     bool    1    0" << std::endl;
              This uses model_MS_Global_a1n_etaa3_* which assumes a1(n,l) = a1(n)
        
        Ex 3: average_a1nl     bool    0    1" << std::endl;
              This uses model_MS_Global_a1l_etaa3_* which assumes a1(n,l) = a1(l)
        
        Ex 4: average_a1nl     bool    0    0
              This uses model_MS_Global_a1nl_etaa3_* which assumes a1(n,l) are all free

### v1.2.2 Minor Changes ### 
	* Corrected issues:
		* Problem: Initial values for labels and units (empty vector of strings) was not allowing more than two columns input data. This is a problem when the user wish to
		            use a 'best_models_fit.ascii' file generated by the IDL_PostMCMC program
        * Solution: The program was modified such that I initialise the labels and units vector with 5 columns instead of 2.

### v1.2.1: Major Changes ### 
    * New functionalities:
   		* The program now uses openmp to parallelise the parallel chains computation. The number of used thread is therefore 
   		   controlled by the OMP_NUM_THREADS=X with X the number of used cpus. Gain are typically optimal when OMP_NUM_THREADS=Nchains/2
   		   Gains are null if OMP_NUM_THREADS>Nchains.
   		* Added scripts for running the program on NYUAD-DALMA (slurm scripting)
   		* Multiples instances can now run using the same executable (required by parallelisation)
   		* The program is now executed as a command line in order to process a specific range of object of the 'ids' list in config_presets.cfg
   		   The command line is of the following syntax: ./tamcmc execute [0/1] [first line number to process] [last line number to process]
   		   Warning : The line id index begins at 0 while the the argument of the command line are begining to 1.
   		   Warning2: The variables 'start_index_ids' and 'last_index_ids' in config_presets.cfg are now obselete
   		             These will be removed on a future version 
   	* Corrected issues:
   		* The clock was not working properly due to parallelisation. This has been rewritten by Benoit Marchand
   		* There was a bug on the models.cpp functions that was fixing the l=0 asymetry to eta (thus to almost 0)
   		* The variable do_backup_cfg_files was creating a backup of the data file while this is controlled by do_backup_input_file
   		
   	* Identified issues:
   		 * Problem: In some rare occasions, the following error appears:
			terminate called after throwing an instance of 'std::ios_base::failure[abi:cxx11]'
			what(): cannot open pipe gunplot -persist: isotream error
		    Reason: The problem is due to memory limits. If the system runs out of memory, gnuplot fails
		    Solution: Reduce Nbuffer so that memory usage remains within your system specifications. A future
		    	      will be more optimised in terms of memory usage so that this kind of situation may be more
		    	      unlikely to face
		    	      
		* Problem: In some rare occasions, 'Nan' are return in the vector of parameters. These models are rejected automatically and do not prevent the proper execution of the code

		* Problem: In some rare occasions, build_lorentzian encouter 'Nan' for some parameters, which make the code to unexpectely stop.
		   Solution: TEMPORARY SOLUTION IS to restart the analysis from where it stopped (restore = 3)

		* Problem: When in restore=3 and if Nbuffer is changed relative to the old run, saving of all variables might not always
		            happens when it should.
		   Reason: It is likely due to the fact the logical test to decide when to save is not accounting well from i0
		   Solution: account for i0 when deciding when to save

		* Problem: Verbose missing when writing the evidence
		   Reason: No verbose yet
		   Solution: Implement a verbose
		   
	    * Problem: No evidence written if all diagnostics are OFF
		   Reason: Investigate the issue

### v1.1.3.4: Minor improvements ### 	
	* Added some verbose and security check in build_lorentzian.cpp. The goal is to see if the bug we encounter once every
	  ~10 Million samples (problem 1 in known existing issues) is coming from there.
	  
	* Corrected issues:
		1. Problem: when restore=3, the acceptance rate plots is not correct.
		   Reason: The exchange file that contains the rejection rate is not defining correctly the x-axis
		   Solution: modify outputs::reject_rate() so that it defines an x-axis which account for the samples that were already processed during the last execution (variable outputs::buf_restore.Nsamples_sofar)

	* Known existing issues:
		* Problem: In some rare occasions, the following error appears:
			tamcmc.out: eigen/Eigen/src/Core/MapBase.h:148: Eigen::MapBase<Derived, 0>::MapBase(Eigen::MapBase<Derived, 0>::PointerType, Eigen::MapBase<Derived, 0>::Index, Eigen::MapBase<Derived, 0>::Index) [with Derived = Eigen::Block<const Eigen::Matrix<double, -1, 1>, -1, 1, false>; Eigen::MapBase<Derived, 0>::PointerType = const double*; Eigen::MapBase<Derived, 0>::Index = long int]: Assertion `(dataPtr == 0) || ( nbRows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == nbRows) && nbCols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == nbCols))' failed.
		   Solution: Contrary to what I was thinking initially, this error seems not related to a plot function. But rather to the truncation approximation of the Lorentzian.
		   REASON UNKNOWN. CIRCUNSTANCES: When truncating the Lorentzians?

		* Problem: In some rare occasions, the following error appears:
			terminate called after throwing an instance of 'std::ios_base::failure[abi:cxx11]'
			what(): cannot open pipe gunplot -persist: isotream error
		   Solution: I posted the problem in stack overflow
		
		* Problem:  The variable do_backup_cfg_files seems also to create a backup of the input file (normally controlled by do_backup_input_fil)
		   Reason: Logic condition is incorrect
                   Solution: Put the correct logic condition

		* Problem: Verbose missing when writing the evidence
		   Reason: No verbose yet
		   Solution: Implement a verbose

		* Problem: Missing parallelisation
		   Solution: Wait that Benoit Lemarchand review the code and implement the parallelisation in it. Deadline given by Benoit is 15/01/2017

###  v1.1.3.3: Minor improvements ### 
	* Corrected issues:
		* Problem: Error messages appear in machines which have gunplay installed, but not gnuplot-qt. 
		   Reason: in diagnostics.cpp, we switch from a ps terminal into a qt terminal to reinitialise the terminal
 		   Solution: remove the line set term 'qt' as it is not necessary to reset the terminal type.
		   tested: Requires to be tested on othxeon12
		* Problem: When handling a list of process it is hard to keep track of which have failed and which have finished and which are ongoing
		   Solution: Generate a simple output log file on the following text format:
			[Number of the process]   [Name of the process]     [processing_status]
		        Use the following codes for processing_status: No processing, Pending, Ongoing , Finished   	 

	* Known existing issues:
		* Problem: In some rare occasions, the following error appears:
			tamcmc.out: eigen/Eigen/src/Core/MapBase.h:148: Eigen::MapBase<Derived, 0>::MapBase(Eigen::MapBase<Derived, 0>::PointerType, Eigen::MapBase<Derived, 0>::Index, Eigen::MapBase<Derived, 0>::Index) [with Derived = Eigen::Block<const Eigen::Matrix<double, -1, 1>, -1, 1, false>; Eigen::MapBase<Derived, 0>::PointerType = const double*; Eigen::MapBase<Derived, 0>::Index = long int]: Assertion `(dataPtr == 0) || ( nbRows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == nbRows) && nbCols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == nbCols))' failed.
		   Solution: Contrary to what I was thinking initially, this error seems not related to a plot function. But rather to the truncation approximation of the Lorentzian.
		   REASON UNKNOWN. CIRCUNSTANCES: When truncating the Lorentzians?

		* Problem: In some rare occasions, the following error appears:
			terminate called after throwing an instance of 'std::ios_base::failure[abi:cxx11]'
			what(): cannot open pipe gunplot -persist: isotream error
		   Solution: I posted the problem in stack overflow
		
		* Problem: when restore=3, the acceptance rate plots is not correct.
		   Reason: The exchange file that contains the rejection rate is not defining correctly the x-axis
		   Solution: modify outputs::reject_rate() so that it defines an x-axis which account for the samples that were already processed during the last execution (variable outputs::buf_restore.Nsamples_sofar)

		* Problem:  The variable do_backup_cfg_files seems also to create a backup of the input file (normally controlled by do_backup_input_fil)
		   Reason: Logic condition is incorrect
                   Solution: Put the correct logic condition

		* Problem: Verbose missing when writing the evidence
		   Reason: No verbose yet
		   Solution: Implement a verbose

		5. Problem: Missing parallelisation
		   Solution: Wait that Benoit Lemarchand review the code and implement the parallelisation in it. Deadline given by Benoit is 15/01/2017


### v1.1.3.2: Minor improvements ###
	* Corrected issues:
		- Fixed an issue with the plots of the models (diagnostics.cpp): The smooth coeficient was improperly calculated such that plots of models were not readable in some cases.
		- Fixed an issue with the plots of the likelihood (diagnostics.cpp): When the likelihood was large, a truncation in the exchange file was making the plots awkward.

	* Added functionalities:
		* Handling the outputs is now made easier:
			- In config_presets.cfg, the new variable 'cfg_out_dir=' allows you to set the main directory for all outputs files (diags, diags/pdfs, outputs, restore). 
			  No need anymore to use config_default.cfg for this: Outputs.output_dir, Outputs.restore_dir, Diagnostics.output_dir are bypassed and superseeded by cfg_out_dir
			- Within cfg_out_dir, output directories are automatically created. The following tree is checked and created/completed is required:
				.../[table_ids(i,0)]			
						.../outputs
						.../diags
						     .../pdfs
						.../restore
			   With table_ids specifies the name of the object that is analysed, according to the table at the bottom of config_preset.cfg
		* [BETA TEST] Saving averaged covariance matrix over Nbuffer or rest (depends how many samples remain)... these might provide more reliable covariance matrix out of the Learning process.
		   REQUIRES THOROUGH TESTING BEFORE VALIDATION

		* Added the possibility to backup the input configurations (*.model, *.cfg) and/or the input data (*.data) into the output directory. The backup is made by zipping the files in 'inputs_backup/'
		   This is controled by two new variables in config_default.cfg, do_backup_cfg_files and do_backup_input_file

		* Added the evidence calculation. Outputs are handled by diagnostics.cpp and are save in the diags directory (text file)

		* Added the possibility of showing/not showing the original data (ie, raw data without smooth) in the plots of models. The control variable is 'show_original_data' in config_default.cfg
		   Setting show_original_data is recommended for very noisy data as it improves the rendering of the plots.

	* Known existing issues:
		* In some rare occasions, the following error appears:
			tamcmc.out: eigen/Eigen/src/Core/MapBase.h:148: Eigen::MapBase<Derived, 0>::MapBase(Eigen::MapBase<Derived, 0>::PointerType, Eigen::MapBase<Derived, 0>::Index, Eigen::MapBase<Derived, 0>::Index) [with Derived = Eigen::Block<const Eigen::Matrix<double, -1, 1>, -1, 1, false>; Eigen::MapBase<Derived, 0>::PointerType = const double*; Eigen::MapBase<Derived, 0>::Index = long int]: Assertion `(dataPtr == 0) || ( nbRows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == nbRows) && nbCols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == nbCols))' failed.
		   This error seems to occurs only in 'nohup' mode and might be related to plot functions. 
		   REASON UNKNOWN. CIRCUNSTANCES: Probably only in 'nohup' mode

		* Missing parallelisation

### v1.1.2: Minor improvements ### 
	* Corrected issues:
		* Problem: In the Misc/IDL_PostMCMC/, the index of the pdf for Height/Width/Amplitude begins at 2.
		   Cause: wrong indexation in the MS_GLOBAL_HEIGHT_L, MS_GLOBAL_AMPLITUDES and MS_GLOBAL_WIDTH
		   Solution: (1) put index=index+1 at the end of the loop and (2) replace index+1 by index when formating the number prior generating the filename
		* Problem: Diagnostic plots are not properly shown (edges missing)
		   Cause: Wrong setup for the margnins in diagnostics.cpp
		   Solution: Change the bottom and left margin from 0.02 (or 0.03) to 0.06

	* Added functionalities:
		* Possibility to switch off/on the smoothness condition by adding 'freq_smoothness  bool  [1/0]   [smoothness coefficient]' in the .model file.
		   Example: 
			- Turn ON the smoothness, with a smoothness coeficient of 1.5 microHz/n^2:  freq_smoothness  bool  1   1.5
			- Turn OFF the smoothness:  freq_smoothness  bool  0   1.5
		* The prior Jeffrey_abs() was added to the list of possible prior. Example of use: For the mode asymetry.

	* Known existing issues:
		* Missing functionality calculating the evidence
		* In some rare occasions, the following error appears:
			tamcmc.out: eigen/Eigen/src/Core/MapBase.h:148: Eigen::MapBase<Derived, 0>::MapBase(Eigen::MapBase<Derived, 0>::PointerType, Eigen::MapBase<Derived, 0>::Index, Eigen::MapBase<Derived, 0>::Index) [with Derived = Eigen::Block<const Eigen::Matrix<double, -1, 1>, -1, 1, false>; Eigen::MapBase<Derived, 0>::PointerType = const double*; Eigen::MapBase<Derived, 0>::Index = long int]: Assertion `(dataPtr == 0) || ( nbRows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == nbRows) && nbCols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == nbCols))' failed.
		   This error seems to occurs only in 'nohup' mode and might be related to plot functions. REASON UNKNOWN. CIRCUNSTANCES: Probably only in 'nohup' mode
		* Saving averaged covariance matrix over Nbuffer... these might provide more reliable covariance matrix out of the Learning process
		* Missing parallelisation


### v1.1.1: Fully functionnal release for Burn-in/Learning/Acquire mode ### 
	* Corrected issues: 
		  * Problem: In the Acquire phase, the sampling proceed for few ~100 samples but then the chain get stuck at one position in the parameter space
			 Cause: Seems to come from the Chain Mixing... If We impose that there is no mixing, then the problem disapear. Further examination of the my program
				shows that only the value of the parameters of the coldest chain were returned.
			 Solution: Impose that the parameters (via params) are initialized to the restore value for all parallel chain (changes in model_def.cpp)
		  * Problem: The initial acceptance rate is rather low for high temperature chains
		  	 Cause: The initial covariance matrix was not properly scaled for chains of T>1. 
		  	 Solution: Corrected by changing in MALA::init_proposal from,
		  	 	initial sigma[m]=std::pow(2.38,2)*(0.1*m+1.)/Nvars into, 
		  	 	sigma[m]=std::pow(2.38,2)*std::pow(Tcoefs[m],0.25)/Nvars
		  	 	This scaling ensures that sigma is always properly initialised accordingly to errors given in errors_default.cfg
		  * Problem: When loading from a restore file, huge quantities of "Warning: Hit epsilon1 in p1"
		  	 Cause: Corrolary of problem 1 and 2. 
		  	 Solution: Solving problem 1 and 2 solves this issue.
		  * Problem: The Learning is made by swapping chains. This should be avoided on the preset mode
		  	 Cause: The dN_mixing variable is updated to dN_mixing = Nsamples + 1
		  	 Solution: Hardcoding of dN_mixing = Nsamples + 1 in Learning mode
		  * Problem: In binary mode, the header does not specify the total number of samples saved so far.
			  Cause: Need to add a line in the header with that information
			  Solution: Added a line in the header file, with keyword 'Nsamples_done'

### v1.1: Second release with added functionalities ### 

### v1.0: Initial release with basic working configuration ### 

	* Known existing issues:
		* Problem with the diagnostics plots: These are not properly shown on some postscript pages. A redefinition of the drawing area is required.
		* Missing functionality calculating the evidence
		* No way to switch ON/OFF the smoothness condition on frequencies
		* Missing parallelisation

	* Added functionalities:
		* In main program: Added an indicator of the swaping rate for the parallel chains
		* In Misc: small routines that can read the binary files have been added (BETA VERSION)
