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

    cout << "[reunion] " << endl;
    // A verifier
    for (int i = 0; i < taille1; i++)
        result[k++] = data1[i];

    for (int i = 0; i < taille2; i++)
        result[k++] = data2[i];

    sort(result, result + k);
    cout << "[reunion] New size = " << k << endl;
}

void partitionH(int pivot, int* data, int &taille,
               int* dataInf,int *taille1,
               int* dataSup,int *taille2) {
    int j = 0;
    int k = 0;

    for (int i = 0; i < taille; i++) {
        //cout << "********************* data[i] = " << data[i] << "; pivot = " << pivot << endl;
        if (data[i] < pivot)
            dataInf[k++] = data[i];
        else
            dataSup[j++] = data[i];
    }
    
    // cout << "[partition] taille1 = " << k << "; taille2 = " << j << endl;

    *taille1 = k;
    *taille2 = j;
}

void exchange(int* data, int* taille, int etape) {
    int myPE;
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    // Find the neighbor
    int neighbor = myPE ^ (0x1 << (etape - 1));
    //int neighbor = myPE ^ (int)pow(2, etape);

    cout << "[exchange] " << "myPE : " << myPE << " / size : " << *taille << " / neighbor : " << neighbor << endl;

    // Send the size to the neighbor
    //cout << "exchange send size : Send by " << neighbor << " : " << myPE << " et etape : " << etape << " size : " << taille << endl;
    MPI_Send(taille, 1, MPI_INT, neighbor, 666, MPI_COMM_WORLD);
    if (*taille > 0) {
        // Send the array to the neighbor
        //cout << "exchange send array : Send by " << neighbor << " : " << myPE << " et etape : " << etape << endl;
        MPI_Send(data, *taille, MPI_INT, neighbor, 667, MPI_COMM_WORLD);
    }

    // Receive size from the neighbor
    //cout << "Neighbo : Recv by " << neighbor << endl;
    MPI_Recv(taille, 1, MPI_INT, neighbor, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (*taille > 0) {
        // Receive array from the neighbor
        MPI_Recv(data, *taille, MPI_INT, neighbor, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    cout << "[exchange] New taille = " << *taille << " (myPE = " << myPE << ")"<<endl;
}

void diffusion(int *pivot, int etape) {
    // Send pivot depending on the etape
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    cout << "[diffusion] : " << "myPE : " << myPE << " / etape : " << etape << endl;

    // Root is the process of broadcast
    int root;
    if (etape == log2(p))
        root = 0;
    else
        root = ((myPE >> etape) << etape);
    
    // cout << "[diffusion] The root is : " << root << endl;
    
    // myPE - root = relative index in sub-hypercude of dimension i
    for (int k = 0; k < etape; k++) {
        if ((myPE - root) < (0x1 << k)) {
            //cout << "Send by " << myPE << " : " << myPE + (0x1 << k) << " et etape : " << etape << endl;
            // Send the pivot
            cout << "[diffusion] : I send the pivot : " << *pivot << " to my friend : " << myPE + (0x1 << k) <<  endl;
            MPI_Send(pivot, 1, MPI_INT, myPE + (0x1 << k), 668, MPI_COMM_WORLD);
        }
        else if ((myPE - root) < (0x1 << (k + 1))) {
            // Receive the pivot
            //cout << "Diffusion : Recv by " << myPE - (0x1 << k) <<  " et etape : " << etape << endl;
            MPI_Recv(pivot, 1, MPI_INT, myPE - (0x1 << k), 668, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            cout << "[diffusion] : I " << myPE << " received a new pivot : " << *pivot << " from my friend : " << myPE - (0x1 << k) << endl;
        }
    }
}

void quickSort(int* data, int& taille) {
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
    // Pretri local
    sort(data, data + taille);

    // Dimension
    int d = log2(p);
    // cout << "Dimension(d) : " << d << endl;

    // Tableau de resultat
    int* dataInf = new int[taille];
    int* dataSup = new int[taille];
    int tailleInf = 0;
    int tailleSup = 0;
    int pivot = 0;

    for (int i = d; i > 0; i--) {
        // Choose the pivot
        if (taille != 0) {
            pivot = data[taille / 2];
            cout << "Pivot[" << i << "] = " << pivot << endl;
        }
        //dataInf = new int[taille];
        //dataSup = new int[taille];

        // Call diffusion
        diffusion(&pivot, i);

        // cout << "         " << myPE << "        New pivot = " << pivot << endl;

        // Partition
        partitionH(pivot, data, taille, dataInf, &tailleInf, dataSup, &tailleSup);

        cout << "I am " << myPE << " and tailleInf = " << tailleInf << "; and tailleSup = " << tailleSup << endl;

        // Exchange of arrays
        if (!(myPE & (0x1 << (i - 1))))
            exchange(dataSup, &tailleSup, i);
        else
            exchange(dataInf, &tailleInf, i);

        // Reunion de listes
        //int result = new int[(int *)(taille) + tailleNeighbor];
        // Si le bit i est a 0
        if (!(myPE & (0x1 << i)))
            reunion(dataInf, tailleInf, dataSup, tailleSup, data);
    }
}


void printAll(int* data,int taille) {  
   int p,myPE;
   MPI_Comm_size(MPI_COMM_WORLD,&p);
   MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
   int *recvcounts, *displs, *recvbuf;
   if (myPE == 0) recvcounts = new int[p];
   MPI_Gather(&taille,1,MPI_INT,recvcounts,1,MPI_INT,0,MPI_COMM_WORLD);
   if (myPE == 0) {
      displs = new int[p];
      displs[0] = 0;
      for (int pe=1;pe<p;pe++) displs[pe] = displs[pe-1]+recvcounts[pe-1];
   }
   if (myPE == 0) recvbuf = new int[displs[p-1]+recvcounts[p-1]];
   MPI_Gatherv(data,taille,MPI_INT,recvbuf,recvcounts,displs,MPI_INT,
               0,MPI_COMM_WORLD);
   if (myPE == 0)
      for (int k=0;k<displs[p-1]+recvcounts[p-1];k++) cout << recvbuf[k] << endl;
   if (myPE == 0) delete recvbuf,recvcounts,displs;
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

/* 
void printAll(int* data, int taille) {
   int p, myPE;
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
   int* res;

   if (myPE == 0)
       res = new int[p * taille];

   MPI_Gather(data, taille, MPI_INT, res, taille, MPI_INT, 0, MPI_COMM_WORLD);

   if (myPE == 0)
      for (int k = 0; k < p * taille; k++) cout << res[k] << endl;

   if (myPE == 0) delete res;
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int myPE;
  MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

  int seed = atoi(argv[1]) + myPE;
  srand(seed);

  int tailleLoc = atoi(argv[2]); // tailleLoc vaut n/p

  int* dataLoc = new int[tailleLoc];

  for (int k = 0; k < tailleLoc; k++)
      dataLoc[k] = rand() % 1000;

  printAll(dataLoc, tailleLoc);

  quickSort(dataLoc, tailleLoc);

  printAll(dataLoc, tailleLoc);

  MPI_Finalize();
  return 0;
}
*/
