
// Random Number not consistent solved
// Just add "#pragma omp critical" to avoid race condition
// Another way would be read the AtomData within openmpi block

#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define CUT_OFF 0.28
#define X 1503.02
#define Y 210.62
#define Z 210.62
#define LEN(x)  (sizeof(x) / sizeof((x)[0]))
#define NBIN 500


struct AtomData 
{
	int noOfAtoms;
	float* atomCoords[3];
}AtomDataObj;

void readAtoms()
{	
	FILE *fp = fopen("new0.25", "r");	
	fscanf(fp, "%d", &AtomDataObj.noOfAtoms);
	int i;
	for (i = 0; i < 3; i++) 
	{
		AtomDataObj.atomCoords[i] = (float*)malloc(AtomDataObj.noOfAtoms*sizeof(float));
	}
	for (i = 0; i < AtomDataObj.noOfAtoms; i++) 
	{
		fscanf(fp, "%f%f%f", &AtomDataObj.atomCoords[0][i], &AtomDataObj.atomCoords[1][i], &AtomDataObj.atomCoords[2][i]);
	}
}

void main(int argc,char **argv) 
{
	int nbin, myid, nproc, nthreads, tid, totalCount=0;
        int sum[100]={0};
	readAtoms();
	MPI_Init(&argc,&argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

        srand((unsigned int)(123456));
	#pragma omp parallel private(tid)
	{
            tid = omp_get_thread_num();
            nthreads = omp_get_num_threads();
	    float x, y, z;
//            srand((unsigned int)(123456));
            int i, j;
            for (i = 0; i < NBIN; i++)
            {
                x = ((float)rand() / (float)(RAND_MAX)) * X;
                y = ((float)rand() / (float)(RAND_MAX)) * Y;
                z = ((float)rand() / (float)(RAND_MAX)) * Z;
                #pragma omp critical
		for (j = 0; j < AtomDataObj.noOfAtoms; j++)
                {
                        float xa = AtomDataObj.atomCoords[0][j];
                        float ya = AtomDataObj.atomCoords[1][j];
                        float za = AtomDataObj.atomCoords[2][j];
                        if ((fabs(x - xa) <= CUT_OFF) && (fabs(y - ya) <= CUT_OFF) && (fabs(z - za) <= CUT_OFF))
                        {
                                printf("\nX: %f, Y: %f, Z: %f ATOM: %d, TID: %d\n", x, y, z, j, tid);
                                sum[tid]++;
                                break;
                        }
                 }
           }
        }
        printf("\nThreads: %d", nthreads);
        for (tid = 0; tid < nthreads; tid++)
        {
              printf("\ntid: %d, count: %d", tid, sum[tid]);
              totalCount += sum[tid];
        }	
//	MPI_Allreduce(&count,&totalCount,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  
	
        printf("\nTotal Count = %d\n",totalCount);
	
	MPI_Finalize();
}
