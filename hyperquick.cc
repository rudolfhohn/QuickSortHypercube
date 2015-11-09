#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <math.h>

using namespace std;

void reunion(int* data1, int taille1,
             int* data2, int taille2, int* result) {
    int k = 0;

    cout << "reunion" << endl;
    // A verifier
    for (int i = 0; i < taille1; i++)
        result[k++] = data1[i];

    for (int i = 0; i < taille2; i++)
        result[k++] = data2[i];

    sort(result, result+k);
}

void partitionH(int pivot, int* data, int taille,
               int* dataInf,int& taille1,
               int* dataSup,int& taille2) {
    int j = 0;
    int k = 0;

    cout << "partitionH" << endl;

    for (int i = 0; i < taille; i++) {
        if (data[i] <= pivot)
            dataInf[k++] = data[i];
        else
            dataSup[j++] = data[i];
    }

    taille1 = k;
    taille2 = j;
}

void exchange(int* data, int& taille, int etape) {
    int myPE;
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    cout << "exchange by : " << myPE << endl;

    // Find the neighbor
    int neighbor = myPE ^ (0x1 << etape);

    // Send the size to the neighbor
    MPI_Send(&taille, 1, MPI_INT, neighbor, 666, MPI_COMM_WORLD);
    if (taille > 0)
        // Send the array to the neighbor
        MPI_Send(data, taille, MPI_INT, neighbor, 667, MPI_COMM_WORLD);

    // Receive size from the neighbor
    MPI_Recv(&taille, 1, MPI_INT, neighbor, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (taille > 0)
        // Receive array from the neighbor
        MPI_Recv(data, taille, MPI_INT, neighbor, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

}

void diffusion(int pivot, int etape) {
    // Send pivot depending on the etape
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    cout << "diffusion by : " << myPE << endl;

    // Root is the process of broadcast
    int root;
    if (etape == pivot)
        root = 0;
    else
        root = ((myPE >> etape) << etape);

    // myPE - root = relative index in sub-hypercude of dimension i
    for (int k = 0; k < etape - 1; k--) {
        if ((myPE - root) < (0x1 << k))
            // Send the pivot
            MPI_Send(&pivot, 1, MPI_INT, myPE + (0x1 << k), 668, MPI_COMM_WORLD);
        else if ((myPE - root) < (0x1 << (k + 1)))
            // Receive the pivot
            MPI_Recv(&pivot, 1, MPI_INT, myPE + (0x1 << (k + 1)), 668, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void quickSort(int* data, int& taille) {
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
    // Pretri local
    sort(data, data + taille);

    // Dimension
    int d = log2(taille);

    // Tableau de resultat
    int* dataInf;
    int* dataSup;
    int tailleInf;
    int tailleSup;


    for (int i = d; i > 0; i--) {
        // Choose the pivot
        int pivot = data[taille / (2 * p)];

        // Call diffusion
        diffusion(pivot, i);

        // Partition
        partitionH(pivot, data, taille, dataInf, tailleInf, dataSup, tailleSup);

        // Exchange of arrays
        if (!(myPE & (0x1 << (i - 1))))
            exchange(dataSup, tailleSup, i);
        else
            exchange(dataInf, tailleInf, i);

        // Reunion de listes
        //int result = new int[(int *)(taille) + tailleNeighbor];
        // Si le bit i est a 0
        //if (!(myPE & (0x1 << i)))
        reunion(dataInf, tailleInf, dataSup, tailleSup, data);
    }
}

void printAll(int* data,int taille) {
   int p,myPE;
   MPI_Comm_size(MPI_COMM_WORLD,&p);
   MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
   int* res;

   if (myPE == 0)
       res = new int[p*taille];

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

  for (int k=0;k<tailleLoc;k++)
      dataLoc[k] = rand()%1000;

  quickSort(dataLoc,tailleLoc);

  printAll(dataLoc,tailleLoc);

  MPI_Finalize();
  return 0;
}
