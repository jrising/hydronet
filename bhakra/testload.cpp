#include "BhakraTest.h"

int main(int argc, const char* argv[])
{
  MPI_Init(NULL,NULL); {
	int rank,size;
	MPI_Comm_rank(MCW,&rank);
	MPI_Comm_size(MCW,&size);

	BhakraTest test;

	vector< GeneticData<map<string, Measure> >* > organisms;

	ifstream isavefile;
	isavefile.open("newsaved.gas");
	if (!isavefile.is_open()) {
	  cout << "===== Base Case =====" << endl;
	  map<string, Measure> params = test.getParameters();
	  double success = .555555;

	  GeneticData<map<string, Measure> >* eve = new GeneticData<map<string, Measure> >(params, 1, success);
	  eve->addOffspring(success, 10);

	  organisms.push_back(eve);

          ofstream osavefile;
          osavefile.open("newsaved.gas");
          osavefile << organisms;
          osavefile.close();
	} else {
	  isavefile >> organisms;
	  isavefile.close();
	}
  }
}
