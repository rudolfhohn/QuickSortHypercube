#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;

void reunion(int* data1,int taille1,
             int* data2,int taille2,int* result) {  }

void partition(int pivot, int* data, int taille,
               int* dataInf,int& taille1,
               int* dataSup,int& taille2) {  }

void exchange(int* data,int& taille,int etape) {  }

void diffusion(int pivot,int etape) {  }

void quickSort(int* data,int& taille) {  }

void printAll(int* data,int taille) {  
   int p,myPE;
   MPI_Comm_size(MPI_COMM_WORLD,&p);
   MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
   int* res;
   if (myPE == 0) res = new int[p*taille];
   MPI_Gather(data,taille,MPI_INT,res,taille,MPI_INT,0,MPI_COMM_WORLD);
   if (myPE == 0) 
      for (int k=0;k<p*taille;k++) cout << res[k] << endl;
   if (myPE == 0) delete res;
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int myPE;
  MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
  int seed = atoi(argv[1])+myPE; 
  srand(seed);
  int tailleLoc = atoi(argv[2]); // tailleLoc vaut n/p
  int* dataLoc = new int[tailleLoc];
  for (int k=0;k<tailleLoc;k++) dataLoc[k] = rand()%1000;  
  quickSort(dataLoc,tailleLoc);
  printAll(dataLoc,tailleLoc);
  MPI_Finalize();
  return 0;
}
