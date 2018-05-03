# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

Programs to perform asteroseismic analysis using a MCMC optinmisation algorithm. The main codes are written in C++. Some additional code is required to perform the postprocessing in e.g. python (not provided here). 

* Version 1.3, read the file changelog.md for a full list of the code history.

### How do I get set up? ###

* All setup is made through configuration files of extension .cfg. So far, three files are required to configure the program: config_presets.cfg, config_default.cfg, errors_default.cfg. Those can be found in the subdirectory Config/.
In addition, an analysis setup file is required (.model) and a data file (.data). All those files are ASCII files.

* Dependencies: Eigen, Boost, cmake,  GSL (optional), OpenMP (optional)

* No exhaustive documentation is yet available, but these examples might be self explanatory. A complete documentation is scheduled to be released soon,

* Deployment instructions: 
Use cmake to install the program. For the main directory of the program:
	mkdir build
	cd build
	cmake ..
	make
Some options exists to deactivate functions (eg. GSL and/or OpenMP)

### Contribution guidelines ###

No external contribution is expected. This project is constantly improved, so please contact me if you need to see some worthy functionnality implemented. 

### Who do I talk to? ###

* Owner: Othman Benomar (NYUAD research associate)
* Contact: ob19@nyu.edu
