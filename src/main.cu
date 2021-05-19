#include <stdio.h>
#include "mbrecon.cuh"
#include <unistd.h>
#include <term.h>
#include <iostream>
#include <string>

using namespace std;

void ClearScreen()
  {
  if (!cur_term)
    {
    int result;
    setupterm( NULL, STDOUT_FILENO, &result );
    if (result <= 0) return;
    }

  putp( tigetstr( "clear" ) );
  return;
	}

// open interactive GUI file selection
string GetFileName()
{
	char filename[1024];
	FILE *f = popen("/usr/bin/zenity --file-selection --modal --title=\"Select file\"", "r");
	fgets(filename, 1024, f);
	string fString (filename);
	fString = fString.substr(0, fString.size() - 1);
	return fString;
}

int main(){

	int userSelection = 1;
	int subSelection = 1;
	mbrecon recon;
	transducerProperties transProp;	

	do{
		ClearScreen();
		cout << ">>> Model-based reconstruction <<<" << endl;
		cout << " [1] Select input dataset" << endl;
		cout << " [2] Define transducer settings" << endl;
		cout << " [3] Define fluence settings" << endl;
		cout << " [4] Define reconstruction settings" << endl;
		cout << " [5] Run reconstruction" << endl;
		cout << " [6] Export functions" << endl;
		cout << " [0] Quit" << endl;
		cin >> userSelection;
		
		if (userSelection == 0){
			// do nothing
		}else if (userSelection == 1)
		{
			recon.load_data(GetFileName());
		}else if (userSelection == 2)
		{
			do{
				ClearScreen();
				cout << ">>> Transducer settings <<<" << endl;
				cout << " > current transducer: " << transProp.getName() << endl; 
				cout << " [1] Print properties" << endl;
				cout << " [2] Modify properties" << endl;
				cout << " [3] Define new" << endl;
				cout << " [4] Save to file" << endl;
				cout << " [5] Load properties from file" << endl;
				cout << " [0] Return" << endl;
				cin >> subSelection; 
				
				if (subSelection == 0)
				{
					// do nothing but leave menu
				}else if (subSelection == 1) // print properties of transducer
				{
					transProp.printProperties();
					cout << "Press ENTER to continue" << endl;
					cin.get();
				}else if (subSelection == 2)
				{
					cout << "Not implemented yet" << endl;
				}else if (subSelection == 3)
				{
					transProp.defineNew();
				}else if (subSelection == 4)
				{
					transProp.saveToFile();
				}else if (subSelection == 5)
				{
					transProp.readFromFile(GetFileName());
				}else
				{
					cout << "Invalid option passed" << endl;
				}

			}while(subSelection != 0);
			recon.set_transProp(transProp);
		}else if (userSelection == 3)
		{

		}else if (userSelection == 4)
		{
			// nothing here yet
		}else if (userSelection == 5){
			recon.recon();
		}else if (userSelection == 6){
			do{
				ClearScreen();
				cout << ">>> Export functions <<<" << endl;
				cout << " [1] Export output volume as vtk" << endl;
				cout << " [2] Export output volume as h5" << endl;
				cout << " [3] Export transducer discretization" << endl;
				cout << " [4] Export transducer sensitivity field" << endl;
				cout << " [0] Back to main menu" << endl;
				cin >> subSelection;

				if (subSelection == 0)
				{
					// do nothing
				}else if (subSelection == 1){
					// TODO check if already reconstructed
					if (recon.get_isRecon())
						recon.export_vtk();
					else
						printf("Cannot export before reconstruction\n");
				}else if (subSelection == 2){
					cout << "Not implemented yet" << endl;
				}
				else{
					cout << "Invlaid selection" << endl;
				}
			}while(subSelection != 0);

		}else{
			cout << "Invalid selection!" << endl;
		}
		

	}while (userSelection != 0);

	return 0;
}
