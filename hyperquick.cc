#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <math>

using namespace std;

void reunion(int* data1,int taille1,
             int* data2,int taille2,int* result) { 
    int k = 0;

    // A verifier
    for (int i = 0; i < taille1; i++) result[k++] = data1[i];
    for (int i = 0; i < taille2; i++) result[k++] = data2[i];

    sort(result, result+k);    
}

void partition(int pivot, int* data, int taille,
               int* dataInf,int& taille1,
               int* dataSup,int& taille2) {
    int j = 0;
    int k = 0;

    for (int i = 0; i < taille; i++) {
        if (data[i] <= pivot)
            dataInf[k++] = data[i];
        else
            dataSup[j++] = data[i];
    }
    
    taille1 = k;
    taille2 = j;
}

void exchange(int* data,int& taille,int etape) { 
     
}

void diffusion(int pivot,int etape) {
    // Il diffuse le pivot selon l'etape ou il est en broadcast
}

void quickSort(int* data,int& taille) { 
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
    // Pretri local
    sort(data, data + taille);

    // Dimension
    int d = log2(taille);

    // Tableau de resultat
    

    for (int i = d; i > 0; i--) {
        // Choix du pivot
        int pivot = data[taille / (2 * p)];
        
        // Appelle diffusion
        diffusion(pivot, 0); 

        // Echange de listes
        int tailleNeighbor;
        exchange(data, &tailleNeighbor, i);
        
        // Reunion de listes
        int result = new int[taille + tailleNeighbor];
        // Si le bit i est a 0
        if (!(myPE & (0x1 << i)))
            reunion(data, taille, data2, taille2, result);
    }
}

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
