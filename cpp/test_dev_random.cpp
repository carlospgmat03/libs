#include <iostream>
using namespace std;
#include <dev_random.cpp>
int main(int argc, char* argv[]) {
	Random semilla;
	if(argc == 1){
		cout << "Random integer using urandom "<<semilla.strong() <<endl ; 
		cout << "Random integer using random  "<< semilla.secure() <<endl<<endl ; 
		cout << "Para generar solo un unsign int ejecute el programa" <<endl;
		cout << "seguido de random: ./test_dev_random random"<<endl;
		cout << endl;
		cout << "Para generar solo un unsign int con urandom ejecute" <<endl;
		cout << "seguido de urandom: ./test_dev_random urandom"<<endl;
		
	}
	else{
		string option=argv[1];
		if(option=="urandom") cout << semilla.strong();
		if(option=="random") cout << semilla.secure();

	}
  return 0;
}
