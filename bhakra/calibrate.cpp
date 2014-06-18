#include "BhakraTest.h"

int main(int argc, const char* argv[])
{
  MPI_Init(NULL,NULL); {
	int rank,size;
	MPI_Comm_rank(MCW,&rank);
	MPI_Comm_size(MCW,&size);

    srand(time(NULL));

    BhakraTest test;
    //test.testSaving("saved.gas", 5);
    test.calibrate(10000, "saved.gas");

    map<string, Measure> params = test.getParameters();
    map<string, Measure>::iterator it;
    for (it = params.begin(); it != params.end(); it++)
      printf("\t%s = %f\n", it->first.c_str(), it->second.getValue());

  } MPI_Finalize();
}
