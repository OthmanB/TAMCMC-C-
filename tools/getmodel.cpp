/*
   A small function that extract one parameter from the binary output
   That parameter is save into a output file
*/

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <iomanip>
#include "models.h"
#include "model_def.h"
#include "data.h"
#include "string_handler.h"
#include "config.h"
#include "version.h"

void showversion();
int options(int argc, char* argv[], int Nmaxlines);
void usage(int argc, char* argv[], int Nmaxlines);
int check_retrocompatibility(VectorXi plength, std::string modelname);
VectorXd adapt2new_MSGlobal(const VectorXi plength, const VectorXd params, const double c0);

Data Data_Nd2Data(Data_Nd dat);

int main(int argc, char* argv[]){

		bool verbose=0;
		int i, Nmaxlines, Nmodels, modelname_switch, testval, status_compat;
		VectorXi plength;
		VectorXd arr0, model0;
		MatrixXd data_out;
		MatrixXd models;
		Data_Nd data;
		Data data_model;
		Config cfg;
		std::string filename_data, filename_params, fileout, modelname;
		std::string char0, line0;
		std::ifstream cfg_session;
		std::ofstream fileout_stream;
		
		Model_def model_list;

		Nmaxlines=5; // Maximum of models that can be provided (and therefore generated)

		testval=options(argc, argv, Nmaxlines);

		filename_data=argv[1];
		filename_params=argv[2];
		modelname=argv[3]; 
	
		if(testval == 10){
			std::cout << " Output file name not provided. Using the default output filename: output_model.ascii" << std::endl;
			fileout="output_model.ascii";
		} else{
			fileout=argv[4];
		}
					
		std::cout << "  0. Configuration: " << std::endl;
		std::cout << "      - Model name: " << modelname << std::endl;		
		std::cout << "      - Data file: " << filename_data << std::endl;
		std::cout << "      - Parameters file: " << filename_params << std::endl;
		std::cout << "      - Output file: " << fileout << std::endl;
	
		std::cout << "  1. Reading the data file..." << std::endl;
		verbose=0;
		data=cfg.read_data_ascii_Ncols(filename_data, " \t", verbose); // the data in a matricial form				
		data_model=Data_Nd2Data(data);	
		modelname_switch=cfg.convert_model_fct_name_to_switch(modelname); // look for the case number that is going to be used in call_model

		std::cout << "  2. Reading the file with the parameters of the model and computing model(s)..." << std::endl;
		cfg_session.open(filename_params.c_str());
		if (cfg_session.is_open()) {
			char0="#";
			std::getline(cfg_session, line0);	
			while(!cfg_session.eof() && char0 == "#"){ // Jump comments lines in the header
				line0=strtrim(line0);
				char0=strtrim(line0.substr(0, 1));
				if (char0 == "#"){
					std::getline(cfg_session, line0);
				}
			}
			// After all the comments, the first line must contain plength
			plength=str_to_Xiarr(line0, " \t"); // The separator is either a space of a tabulation
			status_compat=check_retrocompatibility(plength, modelname); // Verify the consistency of the vector with earlier version of the code

			// All the following lines must correspond to a vector of inputs. 1 line <==> 1 model. Maximum of Nmaxlines lines permitted
			i=0;
			std::getline(cfg_session, line0);
			while(!cfg_session.eof() && i < Nmaxlines){ 
				line0=strtrim(line0);
				char0=strtrim(line0.substr(0, 1));
				
				arr0=str_to_Xdarr(line0, " \t"); // Read the input line
				if(status_compat == 1){ // Attempting adjustement of the size of the vectors
					std::cout << " >>> Incompatibility in vector sizes is attempted to be solved by adding:" << std::endl;
					std::cout << "		(1) c0=20 at the end of the parameter vector" << std::endl;
					std::cout << "		(2) conversion of inclination / splitting into nus.sin(i), nus.cos(i)" << std::endl;
					std::cout << "     This because since version 1.3.0 the Lorentzian truncation parameter is an adjustable variable and the splitting projecions are fitted" << std::endl;
					plength.conservativeResize(plength.size() +1);
					plength[plength.size()-1]=1; // adding c0
					arr0=adapt2new_MSGlobal(plength, arr0, 20.); // last value is c0
				}				
				model0=model_list.call_model_explicit(&data_model, plength, arr0, modelname_switch);

				if(i==0){// Initialisation of the Matrix of parameters
					models.resize(Nmaxlines, model0.size());
				} 
				models.row(i) = model0; //call_model_solo(&data, arr0, plength, modelname);
				std::getline(cfg_session, line0);
				i=i+1;
			}	
			Nmodels=i;	
			cfg_session.close();
 	 	} else {
   			std::cout << "Unable to open the file: " << filename_params << std::endl;
   			std::cout << "Check that the file exist and that the path is correct" << std::endl;
   			std::cout << "Cannot proceed" << std::endl;
   			std::cout << "The program will exit now" << std::endl;
   			exit(EXIT_FAILURE);
		}
		
		std::cout << "  4. Writing outputs into the output file in a simple matricial format... " << std::endl;
		data_out.resize(data.data.rows(), data.data.cols() + Nmodels) ; //We take the initial data and the models as additional columns
		data_out.block(0, 0, data.data.rows(), data.data.cols())=data.data; // Put the data first
		for( i=0; i<Nmodels; i++){
			data_out.col(data.data.cols() + i)=models.row(i);
		}
    	fileout_stream.open(fileout.c_str());
    	if(fileout_stream.is_open()){
			fileout_stream << std::setprecision(12) << data_out << std::endl;
    	} else{
			std::cout << " Unable to open the binary data file " << fileout.c_str() << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
    	}
    	fileout_stream.close();
	
	std::cout << "Output model file successfully written" << std::endl;
}

Data Data_Nd2Data(Data_Nd dat){

		Data dat_out;
		dat_out.x=dat.data.col(0);
		dat_out.y=dat.data.col(1);
		if (dat.data.cols() == 3){
			dat_out.sigma_y=dat.data.col(2);
		}
		dat_out.Nx=dat.data.rows();
	return dat_out;
}

void showversion()
{
    std::cout << APP_NAME " - getmodel tool - " APP_VERSION "\n built on " __DATE__ << std::endl;

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

int options(int argc, char* argv[], int Nmaxlines){

	std::string arg1, arg2;
	int val;
	
	val=-1;
	arg1="";
	
	if(argc == 2){
		arg1=argv[1];
		if(arg1 == "version"){
			 val=0;
		} 
	}
	if(argc == 4){
		val=10;
	}
	if(argc == 5){
		val=11;
	}

	if (val == -1){ // Error code
		usage(argc, argv, Nmaxlines);
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

void usage(int argc, char* argv[], int Nmaxlines){

		if( argc < 4){
			std::cout << " You need to provide at least three argument to that function. The available arguments are: " << std::endl;
			std::cout << "     [1] The data filename. Should be in the same format as those provided to TAMCMC (*.data file)" << std::endl;
			std::cout << "     [2] The filename for the parameters to read. After comments ('#'), this files contains on row(1) plength and on row(2:2+Nmaxlines) the model parameters. Maximum number of models is Nmaxlines=" << Nmaxlines << std::endl;
			std::cout << "         THIS FILE IS NOT A DIRECT OUTPUT OF THE TAMCMC program. You must make it yourself, using e.g. the mean of the parameters (that can be extracted the bin2txt program)" <<std::endl;
			std::cout << "     [3] The model name among the family of MS_Global models (e.g. model_MS_Global_a1etaa3_HarveyLike)" << std::endl;
			std::cout << "     [4] [Optional] the output filename. If not given, then the output file 'output_model.ascii'" << std::endl;
			std::cout << " Call sequence: " << std::endl;
			std::cout << "     " << argv[0] << " [data filename] [model name] [parameter filename] [output filename]" << std::endl;
			exit(EXIT_FAILURE);
		}
	
}

// Small function that verify that plength is of the correct length knowing the modelname. See the full list of referenced modelnames (as listed in config.cpp)
// This function needs a thorough check when problems are detected... 
int check_retrocompatibility(VectorXi plength, std::string modelname){

	int Nplength_expected=-1; // By defaut, consider that plength is not compatible with the current program
	int status=-1;
	if (modelname == "model_MS_Global_a1etaa3_HarveyLike" || modelname == "model_MS_Global_a1etaa3_Harvey1985" || modelname == "model_MS_Global_a1l_etaa3_HarveyLike" ||
	   modelname == "model_MS_Global_a1n_etaa3_HarveyLike" || modelname == "model_MS_Global_a1nl_etaa3_HarveyLike"){
		Nplength_expected = 11;
		if(Nplength_expected != plength.size()){
			std::cout << "   >> Structure of the parameters incorrect. Incompatibility detected" << std::endl;
			std::cout << "   >> Detected size for plength: "<< plength.size() << " ...Expected size:"<< Nplength_expected << " for model '" << modelname << "'" << std::endl;
			if(plength.size() == 10){ // Case where we might know how to get back to our feet
				status=1; 
			} else{
				status=-1; // There is problem that requires a carefull check/debuging
			}
		} else{
			status=2;
		}
	}
	if (modelname == "model_Evolved_Global_a1etaa3_l1mixed"){
		status=2;
	}
	if(status >1){
		std::cout << "   >> Compatibility/Consistency test with earlier version passed..." << std::endl;
	}
	if(status == -1){
		std::cout << "   >> Check \n (1) The compatibility / consistency of your model names with the main program. \n (2) The version of this code matches the version of the main TAMCMC code" << std::endl;
		exit(EXIT_FAILURE);
	}

return status;
}

// Perform the vector adjustements required in case the MS_Global model use the standard predating version 1.3.0
VectorXd adapt2new_MSGlobal(const VectorXi plength, const VectorXd params, const double c0){

	const long double pi = 3.141592653589793238L;
	const int posa1=plength[0]+plength[1]+plength[2]+plength[3]+plength[4]+plength[5]; // old variable a1
	const int posinc=plength.size()-1; // old variable inclination
	const int posa1cosi=plength[0]+plength[1]+plength[2]+plength[3]+plength[4]+plength[5] + 3; // Position of the new variable that are used instead of a1-inclination
	const int posa1sini=plength[0]+plength[1]+plength[2]+plength[3]+plength[4]+plength[5] + 4; // Position of the new variable that are used instead of a1-inclination

	double inclination, a1;	
	VectorXd newparams=params;

	inclination=params[posinc]; // Last elements of the vector was the inclination in version earlier than 1.3.0
	a1=params[posa1]; ; //splitting a1

	newparams[posa1cosi]=sqrt(a1)*cos(inclination*pi/180.);
	newparams[posa1sini]=sqrt(a1)*sin(inclination*pi/180.);

	newparams.conservativeResize(params.size()+1);
	newparams[newparams.size()-1]=c0; //Put a default c0 value at the end of the new vector of parameters

	return newparams;
}
