#include <iostream>
#include <cstring>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#define PI 3.14159265
using namespace std;

//Makes a file "18rings.txt" that represents an initial condition of 18 concentric circles

int main(int argc, char* argv[]) {

	ofstream ofs;
	char test[100] = "";
	char cores[100] = "";
	char p[100] = ""; 
	char labname[100] = ""; 
	char hours[100] = "";
	char mins[100] = "";
	char nodes[100] = "";
	char mem[100] = "";
	vector<string> parameter_names;
	vector<string> parameter_values;
	vector<string> flag_names;
	char DPM[100] = "";
	char TC[100] = "";
	char WUS_CF[100]= "";
	char CK_CF[100] = "";
	char WUS_LEVEL[100] = "";
	char ALT_NAME[100] = "";
	//Include trailing / in string
	string bigdata_path = "/bigdata/alberlab/shared/";
	string subfolder;
	string final_path;
	int divDataCutoff;
	bool bigdata = false;
	bool lab_set = false;
	bool has_subfolder = false;
	bool has_alt_name = false;
	bool printing_vtks = true;
	bool printing_nematic = false;

	for (int i = 1; i < argc; i++) { 
		if (!strcmp(argv[i], "-p")) { 
			strcpy(p,argv[i+1]);
		} else if (!strcmp(argv[i], "-test")) { 
			strcpy(test,argv[i+1]);
		} else if (!strcmp(argv[i], "-hours")) { 
			strcpy(hours,argv[i+1]);
		} else if (!strcmp(argv[i], "-cores")) { 
			strcpy(cores,argv[i+1]);
		} else if (!strcmp(argv[i], "-mins")) { 
			strcpy(mins,argv[i+1]);
		} else if (!strcmp(argv[i], "-nodes")) { 
			strcpy(nodes, argv[i+1]);
		} else if (!strcmp(argv[i], "-mem")) { 
			strcpy(mem, argv[i+1]);
		} else if (!strcmp(argv[i], "-lab")) { 
			lab_set = true;
			strcpy(labname, argv[i+1]);
		} else if (!strcmp(argv[i], "-PRINT")) { 
			printing_vtks = (strcmp(argv[i+1], "0") == 0) ? false : true;
		} else if (!strcmp(argv[i], "-PRINT_NEMATIC")) { 
			printing_vtks = (strcmp(argv[i+1], "1") == 0) ? true : false;
		} else if (!strcmp(argv[i], "-help")) { 
			goto helplabel;
		} else if (!strcmp(argv[i], "-par")) {
			parameter_names.push_back(argv[i+1]);
			parameter_values.push_back(argv[i+2]);
		} else if (!strcmp(argv[i], "-flag")) {
			flag_names.push_back(argv[i+1]);
		} else if (!strcmp(argv[i], "-bigdata")) { 
			bigdata = true;
		} else if (!strcmp(argv[i], "-subfolder")) { 
			has_subfolder = true;
			subfolder = argv[i+1];
		} else if (!strcmp(argv[i], "-sh")) { 
			has_alt_name = true;
			strcpy(ALT_NAME,argv[i+1]);
		} 
	}
	char NAME[100];
	if (has_alt_name) { 
		strcpy(NAME,strcat(ALT_NAME,".sh"));
	} else { 
		strcpy(NAME,"AUTO_BATCH.sh");
	}
	ofs.open(NAME);

	if (!strlen(test)) {
		cout << "Test name: ";
		cin >> test; 
		if (cin.fail()) { 
			cout << "Integer only." << endl;
			ofs.close();
			return 1;
		}
		cin.clear();
	}
	if (!strlen(cores)) { 
		cout << "Number of cores: ";
		cin >> cores; 
		if (cin.fail()) { 
			cout << "Integer only." << endl;
			ofs.close();
			return 1;
		}
	}
	if (!strlen(hours) && !strlen(mins)) { 
		cout << "Hours: ##";
		cin >> hours; 
		if (cin.fail()) { 
			ofs.close();
			return 1;
		}
		cout << "Minutes: ##";
		cin >> mins; 
		if (cin.fail()) { 
			ofs.close();
			return 1;
		}
	} else {
		if (!strlen(hours)) strcpy(hours, "00");
		if (!strlen(mins)) strcpy(mins, "00");
	}
	if (!strlen(p)) {
		cout << "Partition: ";
		cin >> p; 
		if (cin.fail()) { 
			ofs.close();
			return 1;
		}
	}
	if (!strlen(mem)) strcpy(mem, "2");
	if (!strlen(nodes)) strcpy(nodes,"1");


	final_path = (bigdata) ? bigdata_path : "";

	if (has_subfolder) { 
		final_path = final_path + subfolder;
	}


	ofs << "#!/bin/bash -l\n";
	ofs << "#SBATCH --nodes=" << nodes << "\n";
	ofs << "#SBATCH --ntasks=1\n";
	ofs << "#SBATCH --cpus-per-task=" << cores << "\n";
	ofs << "#SBATCH --mem-per-cpu=" << mem << "G\n";
	ofs << "#SBATCH --time=0-" << hours << ":" << mins << ":00\n";
	ofs << "#SBATCH --output=my"<< test << ".stdout\n";
	ofs << "#SBATCH --job-name=\"" << test << "\"\n";
	ofs << "#SBATCH -p " << p << " \n";
	if (lab_set) { 
		ofs << "#SBATCH -A " << labname << " \n";
	}
	ofs << "#SBATCH -p " << p << " \n";
	//ofs << "SBATCH whatever thing for reservation... 
	ofs << "export OMP_NUM_THREADS=" << cores << "\n";
	if (printing_vtks) ofs << "mkdir " << final_path << "Animate_Cyt_" << test << "\n";
	if (printing_nematic) ofs << "mkdir " << final_path << "Nematic_test_" << test <<"\n";
	ofs << "mkdir " << final_path << "Locations_test_" << test << "\n";
	if (printing_vtks) ofs << "mkdir " << final_path << "Animate_No_Cyt_" << test << "\n";
	ofs << "mkdir " << final_path << "Cell_Data_" << test << "\n";
	ofs << "mkdir " << final_path << "Tissue_Data_" << test << "\n";
	ofs << "./program " << final_path << "Animate_Cyt_" << test << " " 
		<< final_path << "Locations_test_" << test << " " 
		<< final_path << "Nematic_test_" << test << " " 
		<< final_path << "Animate_No_Cyt_" << test << " "
		<< final_path << "Cell_Data_" << test << " "
		<< final_path << "Tissue_Data_" << test;
	for (unsigned int i = 0; i < parameter_values.size(); i++ ) { 
		ofs << " " << parameter_names.at(i) << " " << parameter_values.at(i);
	}
	for (unsigned int i = 0; i < flag_names.size(); i++ ) { 
		ofs << " " << flag_names.at(i);
	}
	ofs << endl;
	ofs.close();


	return 0;


	helplabel:
	cout << "This file generates a SBATCH script called AUTO_BATCH.sh.\n";
	cout << "Options and their defaults are:" << endl;
	cout << "-p <partition>.  Omitting this will give a console prompt." << endl;
	cout << "-test <name>.  Ommitting this will give a console prompt." << endl;
	cout << "-hours <##> and -mins <##>.  Two digits only.  Filling only of these" << 
		" will default the other to 00." << endl;
	cout << "-cores <#>. Ommitting this will give a console prompt to define it." <<endl;
	cout << "-nodes <#>. Defaults to 1." << endl;
	cout << "-mem <#>.  Specifies number of gigabytes of memory per CPU. Defaults to 2." << endl;
	cout << "Example:  ./a.out -p short -hours 02 -cores 12 -test mech_div_5" << endl;

	return 0;
}
