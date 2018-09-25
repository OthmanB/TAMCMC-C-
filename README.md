# README #

---
**NOTE FOR PLATO WP 128**

The master version does not support mixed modes fitting. Please use the dev-mixedmodes version instead. 
---

### TAMCMC C++ repository ###

Programs to perform asteroseismic analysis using a MCMC optinmisation algorithm. The main codes are written in C++. Some additional code is required to perform the postprocessing in e.g. python (not provided here). 
The algorithm implementation is based on the IDL implementation of my code described in Benomar et al. (2009), see http://adsabs.harvard.edu/abs/2009A%26A...506...15B
This implementation is itself based on the Metropolis-Hasting-Langevin scheme described in Atchade Y. (2006), see https://link.springer.com/article/10.1007/s11009-006-8550-0

* Version 1.3.1, read the file CHANGELOG.md for a full list of the code history.
	The current implementation focuses on the asteroseismic analysis (individual pulsations) of main sequence stars. For other stars, one would need to write
	dedicated functions, but the code is modular so that only minimal effort is required to do so. 
	
* Future implementations should incorporate C++ version of my codes (currently in Python or IDL) to prepare the configuration files (*.model and *.data)

### How do I get set up? ###

An exhaustive help is available on the wiki of this project (https://github.com/OthmanB/TAMCMC-C/wiki). This readme provides basic explanations only

* All setup is made through configuration files of extension .cfg. So far, three files are required to configure the program: config_presets.cfg, config_default.cfg, errors_default.cfg. Those can be found in the subdirectory Config/.
In addition, an analysis setup file is required (.model) and a data file (.data). All those files are ASCII files.

* Dependencies: Eigen, Boost, cmake,  GSL (optional), OpenMP (optional)


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
