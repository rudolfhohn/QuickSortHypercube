#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <unistd.h>
#include <math.h>

using namespace std;


/**
 * @brief Gather values from two arrays and put them into another one (sorted)
 *
 * @param data1     first array
 * @param taille1   size of first array
 * @param data2     second array
 * @param taille2   size of second array
 * @param result    array to store all values
 * @param taille    size of result (updated)
 */
void reunion(int* data1, int taille1,
             int* data2, int taille2, int*& result, int& taille) {
    int myPE;
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    int k = 0;
    // Update the size
    taille = taille1 + taille2;

    // Clear the array and reallocate the memory with the right size
    delete[] result;
    result = new int[taille];

    // Get values from data1
    for (int i = 0; i < taille1; i++)
        result[k++] = data1[i];

    // Get values from data2
    for (int i = 0; i < taille2; i++)
        result[k++] = data2[i];

    cout << "[reunion " << myPE << "] k = " << (taille1+taille2) << endl;
    // Sort data we just gather
    sort(result, result + taille1 + taille2);
}


/**
 * @brief Depending on the pivot, dispatch value on data inf / sup
 *
 * @param pivot     pivot, value smaller will be in dataInf and higher dataSup
 * @param data      data with all the values
 * @param taille    size of data
 * @param dataInf   array to contain smaller values than the pivot
 * @param taille1   size of dataInf
 * @param dataSup   array to contain higher values than the pivot
 * @param taille2   size of dataSup
 */
void partitionH(int pivot, int* data, int taille,
               int*& dataInf, int& taille1,
               int*& dataSup, int& taille2) {
    int myPE;
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
    int j = 0;
    int k = 0;

    // Find out the size of arrays (inf / sup)
    for (int i = 0; i < taille; i++) {
        if (data[i] < pivot)
            k++;
        else
            j++;
    }

    // Clear arrays and reallocate with the right size
    delete[] dataInf;
    delete[] dataSup;
    dataInf = new int[k];
    dataSup = new int[j];

    k = 0;
    j = 0;
    // Puts values in arrays
    for (int i = 0; i < taille; i++) {
        if (data[i] < pivot)
            dataInf[k++] = data[i];
        else
            dataSup[j++] = data[i];
    }

    // Update the sizes
    taille1 = k;
    taille2 = j;
}


/**
 * @brief Function to exchange two arrays
 *
 * Calculate the neighbor,
 * then send the size and the data.
 *
 * Receive the size from the neighbor and then the data.
 * With the size, we can reallocate the correct size for the data array.
 *
 * @param data      data to exchange
 * @param taille    size of data
 * @param etape     etape of the algo
 */
void exchange(int*& data, int& taille, int etape) {
    int myPE;
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    // Find the neighbor
    int neighbor = myPE ^ (0x1 << etape - 1);

    cout << "[exchange] " << "myPE : " << myPE << " / size : " << taille << " / neighbor : " << neighbor << endl;

    // Send the size to the neighbor
    MPI_Send(&taille, 1, MPI_INT, neighbor, 666, MPI_COMM_WORLD);
    if (taille > 0) {
        // Send the array to the neighbor
        MPI_Send(data, taille, MPI_INT, neighbor, 667, MPI_COMM_WORLD);
    }


    // Receive size from the neighbor
    MPI_Recv(&taille, 1, MPI_INT, neighbor, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (taille > 0) {
        // Clear the data and set the new size
        delete[] data;
        data = new int[taille];
        // Receive array from the neighbor
        MPI_Recv(data, taille, MPI_INT, neighbor, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

}


/**
 * @brief TODO
 *
 * @param pivot Pivot to send
 * @param etape Etape of the algo
 */
void diffusion(int& pivot, int etape) {
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

    // myPE - root = relative index in sub-hypercube of dimension i
    for (int k = 0; k < etape; k++) {
        if ((myPE - root) < (0x1 << k)) {
            cout << "[diffusion] : I send the pivot : " << pivot << " to my friend : " << myPE + (0x1 << k) << endl;
            // Send the pivot
            MPI_Send(&pivot, 1, MPI_INT, myPE + (0x1 << k), 668, MPI_COMM_WORLD);
        }
        else if ((myPE - root) < (0x1 << (k + 1))) {
            // Receive the pivot
            MPI_Recv(&pivot, 1, MPI_INT, myPE - (0x1 << k), 668, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            cout << "[diffusion] : I " << myPE << " received a new pivot : " << pivot << " from my friend : " << myPE - (0x1 << k) << endl;
        }
    }
}


/**
 * @brief Quicksort implementation
 *
 * @param data      data to sort
 * @param taille    size of data
 */
void quickSort(int*& data, int& taille) {
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
    // local sort
    sort(data, data + taille);

    // Dimension
    int d = log2(p);
    cout << "Dimension(d) : " << d << endl;

    // Array with result
    // Memory is not allocated yet,
    // the allocation will be done when we know the size to allocate
    int* dataInf = new int;
    int* dataSup = new int;
    int tailleInf;
    int tailleSup;
    int pivot = 0;

    for (int i = d; i > 0; i--) {
        // Choose the pivot
        if (taille != 0) {
            pivot = data[taille / 2];
            cout << "Pivot[" << i << "] = " << pivot << endl;
        }

        // Call diffusion
        diffusion(pivot, i);

        // Partition
        partitionH(pivot, data, taille, dataInf, tailleInf, dataSup, tailleSup);

        cout << "[Quicksort] PE : " << myPE << " tailleInf : " << tailleInf << endl;
        cout << "[Quicksort] PE : " << myPE << " tailleSup : " << tailleSup << endl;

        // Exchange of arrays
        if (!(myPE >> (i - 1) & 0x1))
            exchange(dataSup, tailleSup, i);
        else
            exchange(dataInf, tailleInf, i);

        reunion(dataInf, tailleInf, dataSup, tailleSup, data, taille);
    }
}


/**
 * @brief Print data
 *
 * @param data      array to print
 * @param taille    size of array
 */
void printAll(int* data,int taille) {
    int p, myPE;

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    int *recvcounts, *displs, *recvbuf;

    if (myPE == 0) recvcounts = new int[p];

    MPI_Gather(&taille, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (myPE == 0) {
       displs = new int[p];
       displs[0] = 0;
       for (int pe = 1; pe < p; pe++) displs[pe] = displs[pe - 1] + recvcounts[pe - 1];
    }

    if (myPE == 0) recvbuf = new int[displs[p - 1] + recvcounts[p - 1]];

    MPI_Gatherv(data, taille, MPI_INT, recvbuf, recvcounts, displs, MPI_INT,
                0,MPI_COMM_WORLD);

    if (myPE == 0)
       for (int k = 0; k < displs[p - 1] + recvcounts[p - 1]; k++) cout << recvbuf[k] << endl;

    if (myPE == 0) delete recvbuf, recvcounts, displs;
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int nbPE, myPE;

    MPI_Comm_size(MPI_COMM_WORLD, &nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    int seed = atoi(argv[1]) + myPE;

    srand(seed);

    int tailleLoc = atoi(argv[2]); // tailleLoc vaut n/p
    int* dataLoc = new int[tailleLoc];

    for (int k = 0; k < tailleLoc; k++) dataLoc[k] = rand() % 1000;

    quickSort(dataLoc, tailleLoc);

    printAll(dataLoc, tailleLoc);

    MPI_Finalize();
    return 0;
}

