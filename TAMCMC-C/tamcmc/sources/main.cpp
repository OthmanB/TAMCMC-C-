/*
 * main.cpp 
 * Created on: 02 Mar 2016
*/

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include "unistd.h"
#include "MALA.h"
#include "config_presets.h"
#include "version.h"

void zip_config(const std::string cfg_file_default, const std::string error_file_default, const std::string cfg_file_presets, 
		const std::string model_file, const std::string output_dir, const std::string id);
void zip_data(const std::string data_file, const std::string output_dir, const std::string id);
std::string shell_exec(const std::string cmd);
bool isdir(const std::string pathname);
bool generate_dir_tree(const std::string rootdir, const std::vector<std::string> subdirs);
void status_logfile( const std::vector< std::vector<std::string> > ids, const int current_proc_ind, const int first_proc_ind, const int last_proc_ind, const std::string logfile);
void showversion();
int  options(int argc, char* argv[]);
void usage(int argc, char* argv[]);

extern long ben_clock();

int main(int argc, char* argv[]){

	bool execute=1; // if argv[1]="execute" then look for argv[2]==execute. 
			// argv[3]=start_index_ids argv[4]=last_index_ids. 
			// If execute = 0 ==> Show only the setup
	int id_index=0, process_index=0,startV=-1,lastV=-1;
	std::string dir_backup;
	int msg_code;
	
	msg_code=options(argc, argv);
	if(msg_code == -1){
		std::cout << "Error detected in options. Cannot proceed. Debug required." << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		if(msg_code == 10){
			std::istringstream(argv[2]) >> execute;
		}
		if(msg_code == 11){
			std::istringstream(argv[2]) >> execute;
			std::istringstream(argv[3]) >> startV;
			std::istringstream(argv[4]) >> lastV;
			startV=startV-1;
			lastV=lastV-1;	
		}
	}

	std::cout << " --------- TAMCMC-CPP ----------" << std::endl;

	// Set the current path variables
	std::string cpath=getcwd(NULL, 0);
	std::cout << "- Current directory: " << cpath << std::endl;

    std::string error_file_default=cpath + "/Config/default/errors_default.cfg"; // This is the file that defines the error values used to initialize the covariance matrix
    std::string cfg_file_default=cpath + "/Config/default/config_default.cfg"; // This is all the parameters required to run the TAMCMC
	
	std::string cfg_file_modelslist=cpath + "/Config/default/models_ctrl.list";
	std::string cfg_file_priorslist=cpath + "/Config/default/priors_ctrl.list";
	std::string cfg_file_likelihoodslist=cpath + "/Config/default/likelihoods_ctrl.list";
	
	std::string cfg_file_presets=cpath + "/Config/config_presets.cfg"; // This is a simpler configuration that allow scripting of the TAMCMC... This is the main configuration file
    std::string logfile=cpath + "/processes.log";
    
	// Load the default configuration
	std::cout << "- Loading the default configurations..." << std::endl;
        Config config(cpath, cfg_file_default, error_file_default, cfg_file_modelslist, cfg_file_priorslist, cfg_file_likelihoodslist);

	// Load the Preset configuration (Master configurator)
	std::cout << "- Loading the Preset configuration: config_presets.cfg..." << std::endl;
	Config_presets config_master(cfg_file_presets, &config); // read the master config file + initialize counters to 0

	if(execute == 0){
		std::cout << "Detection of execute=0... The code will exit now"<< std::endl;
		exit(EXIT_SUCCESS);
	}
    
	long begin_time;
	if(config_master.force_manual_config == 0){
		// BENOIT: after config_preset reset first/last start index to command line parameter values
                config_master.first_id_ind=startV;
                config_master.last_id_ind=lastV;
		// If the user-defined last index is greater than the size of the table, force the last index to match the size of the array
		if (config_master.last_process_ind > config_master.processing.size()){ 
			config_master.last_process_ind=config_master.processing.size()-1;
		}
		// If the user-defined last index is greater than the size of the table, force the last index to match the size of the array
		if (config_master.last_id_ind > config_master.table_ids.size()){
			config_master.last_id_ind=config_master.table_ids.size()-1;
		}
        // ---- Initialize the logfile ----
       // //status_logfile(config_master.table_ids, config_master.first_id_ind, config_master.first_id_ind, config_master.last_id_ind, logfile);

        // ---- Begin to process all the requested objects ----
        for(int i=config_master.first_id_ind; i<=config_master.last_id_ind; i++){ // For each star
            // Updating the logfile after finishing the given object i
           // status_logfile(config_master.table_ids, i, config_master.first_id_ind, config_master.last_id_ind, logfile);
			for(int jj=config_master.first_process_ind; jj<=config_master.last_process_ind; jj++){  // and for each phase
					std::cout << "---------------------------------------------------------------------------------------\n" << std::endl;
					begin_time = ben_clock();
					config_master.current_id_ind=i; // index pointing to the current object id that we have to process
					config_master.current_process_ind=jj;

					std::cout << "                       --------------------------------" << std::endl;
					std::cout << "                       Processing Object " << i+1 << "/" << config_master.table_ids.size() << ": ";
					std::cout << config_master.table_ids[i].at(0) << std::endl;
					std::cout << "                                   Phase " << jj+1 << "/" << config_master.processing.size() << ": ";
					std::cout << config_master.processing[jj] << std::endl;
					std::cout << "                       --------------------------------" << std::endl;

					std::cout << "- Modify the default configuration according to the presets..." << std::endl;
					config_master.apply_presets(&config);

					// Setup the configuration according to user-requested options (either presets configuration or manual configuration)
					std::cout << " Loading the observational constraints (data) and restore setup..." << std::endl;
					config.setup();

					//dir_backup=config.outputs.dir_out + "/" + config_master.table_ids[i].at(0) + "/inputs_backup/";
					if( config.outputs.do_backup_cfg_files == 1){
						std::cout << "    do_backup_cfg_files = 1 ===> Backup of the files *.cfg and *.model in progress..." << std::endl;
						std::cout << "    Target directory: " << dir_backup << std::endl;
						zip_config(cfg_file_default, error_file_default, cfg_file_presets, config.modeling.cfg_model_file, 
							   config_master.cfg_out_dir, config_master.table_ids[i].at(0));
					} else{
						std::cout << "    do_backup_cfg_files = 0 ===> No backup of the files *.cfg and *.model" << std::endl;
					}
				        if( config.outputs.do_backup_input_file == 1){
						std::cout << "    do_backup_input_file = 1 ===> Backup of the files *.data in progress..." << std::endl;
						zip_data(config.data.data_file,  
							 config_master.cfg_out_dir, config_master.table_ids[i].at(0));
					} else{
						std::cout << "    do_backup_input_file = 0 ===> No backup of the file *.data" << std::endl;
					}
					std::cout << std::endl;

					// process the model of the star[i] and for phase[j]... to do so, we need to reinitialize the configuration according to the new presets
					MALA *TAMCMC= new MALA(&config); // Nsamples, Nchains, Nvars
					//std::cout << "Before model_current initialisation" << std::endl;
					Model_def *model_current= new Model_def(&config, TAMCMC->Tcoefs, 1); // Verbose=1 ==> Show variables/constant and consistency test of likelihood/prior 
					//std::cout << "After model_current initialisation" << std::endl;
					Model_def *model_propose= new Model_def(&config, TAMCMC->Tcoefs, 0); // Verbose=0 ==> No need to Show/test anything: already done by model_current
					//std::cout << "After model_propose initialisation" << std::endl;
					Outputs *out = new Outputs(&config, TAMCMC->Tcoefs);
					//std::cout << "After Outputs initialisation" << std::endl;
					Diagnostics *diags = new Diagnostics(&config);

					long begin_time = ben_clock();
					TAMCMC->execute(model_current, model_propose, &config.data.data, &config, out, diags);
					std::cout << " ----------------------------" << std::endl;
					std::cout << "    Calculation finished in: ";
					std::cout << float( ben_clock () - begin_time ) /  60. << " min = ";
					std::cout << float( ben_clock () - begin_time ) /  3600. << " hours = ";
					std::cout << float( ben_clock () - begin_time ) /  86400. << " days" << std::endl;
					std::cout << " ----------------------------" << std::endl;
					std::cout << " ------------------------------------------------------------------------------------" << std::endl;
                
					// Cleaning all the variables to avoid memory dynamic array issues (e.g. memory leaks)
					//config.reset();
					delete TAMCMC;
					delete model_current;
					delete model_propose;
					delete out;
					delete diags;
			}
		}
        // ---- Update one last time the log file ----
       // status_logfile(config_master.table_ids, config_master.last_id_ind+1, config_master.first_id_ind, config_master.last_id_ind, logfile);
	}
	
}

void status_logfile( const std::vector< std::vector<std::string> > ids, const int current_proc_ind, const int first_proc_ind, const int last_proc_ind, const std::string logfile){
    
    std::string status, space;
    std::ostringstream strg;
    std::string filename_log;
    std::ofstream outfile_log;

    filename_log=logfile;
    
    /////// Write the log file ////////
    outfile_log.open(logfile.c_str()); // Overwrite any existing file in PLAIN ASCII
    if (outfile_log.is_open()){
        outfile_log << "# This is a log-file generated by the main program. \n";
        outfile_log << "# This file contains the status for each process from the process list\n" ;
        outfile_log << "# Total number of processes= " << ids.size() << "\n";
        outfile_log << "# Variables names: \n";
        outfile_log << "#      col[0] : index of the process \n";
        outfile_log << "#      col[1] : name of the process \n";
        //outfile_log << "#      col[2] : phase of the process \n";
        outfile_log << "#      col[2] : status for the process. Could be 'No processing', 'Pending', 'Ongoing' or 'Finished' \n";
        
        strg.str(std::string());
        for (int i0=0; i0<ids.size(); i0++){
            status=" ";
            if (i0 >= first_proc_ind && i0 < current_proc_ind){
                status="   Finished     ";
            }
            if (i0 == current_proc_ind){
                status="   Ongoing      ";
            }
            if (i0 > current_proc_ind && i0 <= last_proc_ind){
                status="   Pending      ";
            }
            if (i0 < first_proc_ind || i0 > last_proc_ind){
                status="   No processing";
            }
            if (i0 < 10){
                space="    ";
            }
            if (i0 >= 10 && i0 < 100){
                space="   ";
            }
            if (i0>=100){
                space="  ";
            }
            strg << i0 << space << ids[i0].at(0) << status << "\n";
        }
        outfile_log << strg.str().c_str();
        outfile_log.flush(); // Explicitly specify to flush the data into the disk
        strg.str(std::string()); // clear the strg buffer
        outfile_log.close();
    }
    else {
        std::cout << " Unable to open file " << logfile << std::endl;
        std::cout << " Check that the full path exists" << std::endl;
        std::cout << " The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    }

    
}

void zip_config(const std::string cfg_file_default, const std::string error_file_default, const std::string cfg_file_presets, 
		const std::string model_file, const std::string rootdir, const std::string id){
// 
// This function (1) create a subdirectory inputs_backup into the output_dir if necessary and (2) zip the configuration files into that new directory
//	
	bool checkdir;
	std::string cmd, rcmd, output_dir;
	std::vector<std::string> subdirs;

	subdirs.push_back("/" + id);
	subdirs.push_back("/" + id + "/inputs_backup/");
	generate_dir_tree(rootdir, subdirs);

	output_dir= rootdir + "/" + id + "/inputs_backup/";

/*	checkdir=isdir(output_dir);
	if(checkdir == 0){ // Need to create the directory before exporting the configuration into the zip file
		cmd="mkdir " + output_dir;
		rcmd=shell_exec(cmd);
		std::cout << rcmd << std::endl;
	}
*/
	cmd="zip -9 -j " + output_dir + "cfg_backup.zip " + cfg_file_default + " " + error_file_default + " " + cfg_file_presets + " " + model_file;
	shell_exec(cmd);
	std::cout << rcmd << std::endl;
	std::cout << "Configuration zipped and saved. Check for error messages above..." << std::endl;
}

void zip_data(const std::string data_file, const std::string rootdir, const std::string id){
// 
// This function (1) create a subdirectory inputs_backup into the output_dir if necessary and (2) zip the data files into that new directory
//	
	bool checkdir;
	std::string cmd, rcmd, output_dir;
	std::vector<std::string> subdirs;

	subdirs.push_back(id);
	subdirs.push_back(id + "/inputs_backup/");
	generate_dir_tree(rootdir, subdirs);

	output_dir= rootdir + "/" + id + "/inputs_backup/";

/*	checkdir=isdir(output_dir);
	if(checkdir == 0){ // Need to create the directory before exporting the configuration into the zip file
		cmd="mkdir " + output_dir;
		rcmd=shell_exec(cmd);
		std::cout << rcmd << std::endl;
	}
*/
	output_dir = rootdir + subdirs[1];
	cmd="zip -9 -j " + output_dir + "data_backup.zip " + data_file;
	shell_exec(cmd);
	std::cout << rcmd << std::endl;
	std::cout << "Data file zip and saved. Check for error messages above..." << std::endl;
}


std::string shell_exec(const std::string cmd){
// Small function that execute a script on the shell
// Returns the result of the execution on the form of a string. DUPLICATED IN CONFIG_PRESETS.CPP
    char buffer[128];
    std::string result = "";
    //FILE* pipe = popen(cmd, "r");
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    try {
        while (!feof(pipe)) {
            if (fgets(buffer, 128, pipe) != NULL)
                result += buffer;
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    return result;
}

bool isdir(const std::string pathname){
// Small function that test whether a directory exists
// Return 0 if not, 1 if exists. DUPLICATED IN CONFIG_PRESETS.CPP

	struct stat sb;

	if (stat(pathname.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
	{
	    return 1;
	} else {
	    return 0;
	}
}


bool generate_dir_tree(const std::string rootdir, const std::vector<std::string> subdirs){
// Small function that look into the tree of subdirs whether any directory does not exist
// If it is not, it creates it. DUPLICATED IN CONFIG_PRESETS.CPP

	bool v, vout;
	std::string path, cmd, rcmd;

	vout=1;
	for(int i=0; i<subdirs.size(); i++){
	        path=rootdir + subdirs[i];
		v=isdir(path);
		if(v == 0){
			cmd="mkdir " + path;
			rcmd=shell_exec(cmd);
			std::cout << rcmd << std::endl;
		} else{
			std::cout << path << " already exists " << std::endl;
		}
		vout=v*vout; // If all dirs exist, vout should be 1. Otherwise 0.
	}

return vout;
}


void showversion()
{
    std::cout << APP_NAME " " APP_VERSION "\n built on " __DATE__ << std::endl;

#   if defined(__clang__)
    	printf(" with clang " __clang_version__);
#   elif defined(__GNUC__)
    	printf(" with GCC");
    	printf(" %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#   elif defined(_MSC_VER)
    	printf(" with MSVC");
    	printf(" %d", MSVC_VERSION);
#   else
    printf(" unknown compiler");
#   endif

    std::cout << "\n features:";
#   if defined(__i386__) || defined(_M_IX86)
    std::cout << " i386" << std::endl;
#   elif defined(__x86_64__) || defined(_M_AMD64)
    std::cout << " x86_64" << std::endl;
#   endif
    std::cout << " Author: " << APP_COPYRIGHT << std::endl;

}


int options(int argc, char* argv[]){

	std::string arg1, arg2;
	int val;
	
	val=-1;
	arg1="";
	
	if(argc > 1){
		arg1=argv[1];
	}
	if(argc > 2){
		arg2=argv[2];
	}
	if(argc == 2){
		if(arg1 == "version"){
			 val=0;
		} 
	}
	if(argc == 3){
		if(arg1 == "execute" && arg2 == "0"){
				val=10;
		}
		if(arg1 == "execute" && arg2 == "1"){ // missing 2 arguments
			val=-1;
		}
	}
	if(argc == 5){
		val=11;
	}

	if (val == -1){ // Error code
		usage(argc, argv);
	} 
	if (val == 0){ // Version code
		showversion(); 
		exit(EXIT_SUCCESS);
	}
	if (val > 0 ){
		return val; // Execution code val
	} else{
		return -1; // Default value is to return an error code
	}
}

void usage(int argc, char* argv[]){

			std::cout << "Unrecognized argument" << std::endl;
			std::cout << "     - To execute: " << argv[0] << " execute 1 <start_idx> <last_idx>" << std::endl;
			std::cout << "     - To stop after reading the configuration: " << argv[0] << " execute 0" << std::endl;
			std::cout << "     - To show version: " << argv[0] << " version" << std::endl;
			exit(EXIT_FAILURE);
}
