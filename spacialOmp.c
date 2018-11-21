#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>



void main(int argc,char **argv) 
{
	double t1, t2;
	int nAtoms, myid, nproc, nthreads, tid, totalCount=0, globalCount=0;
	int i, NBIN;
        int sum[100]={0};
	int NX = 5, NY = 2, NZ = 2;
	double X = 1060.28, Y = 212.06, Z = 212.06;
	int GBIN = 7000000;
	double CUT_OFF = 1.45;
	double LX = X/NX;
	double LY = Y/NY;
        double LZ = Z/NZ;	
//	printf("LX: %f, LY: %f, LZ: %f\n", LX, LY, LZ);
        double* atomCoords[3];
        FILE *fp = fopen("50%", "r");
        fscanf(fp, "%d", &nAtoms);
        for (i = 0; i < 3; i++)
        {
                atomCoords[i] = (double*)malloc(nAtoms * sizeof(double));
        }
	for (i = 0; i < nAtoms; i++)
        {
           fscanf(fp, "%lf%lf%lf", &atomCoords[0][i], &atomCoords[1][i], &atomCoords[2][i]);
	}

	MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	NBIN = GBIN / nproc;
	t1 = MPI_Wtime();

	int ii, jj, kk;
        int xn, yn, zn, sid, pointer;
	unsigned int seed;
        double randx, randy, randz;
	int nbin = NBIN / 20;
	double* L_atomCoords[3];

        omp_set_num_threads(20);
	#pragma omp parallel private(tid, xn, yn, zn, sid, pointer, seed, ii, jj, kk, randx, randy, randz, L_atomCoords)
	{
            tid = omp_get_thread_num();
            nthreads = omp_get_num_threads();
	    pointer = 0;

	    for (ii = 0; ii < 3; ii++)
            {
                L_atomCoords[ii] = (double*)malloc(nAtoms * sizeof(double));
            }
	    for (ii = 0; ii < nAtoms; ii++)
       	    {
                 xn = (int) (atomCoords[0][ii] / LX);
		 yn = (int) (atomCoords[1][ii] / LY); 
		 zn = (int) (atomCoords[2][ii] / LZ);
		 sid = xn * NY * NZ + yn * NZ + zn;
		if (sid == tid)
		 {
		    L_atomCoords[0][pointer] = atomCoords[0][ii];
                    L_atomCoords[1][pointer] = atomCoords[1][ii];               
		    L_atomCoords[2][pointer] = atomCoords[2][ii];
		    pointer ++;
		 }
	    }		

            srand((unsigned int)(12345+tid*123+321*myid));
            xn = tid / 4;
            yn = (tid / 2) % 2;
            zn = tid % 2;
	    for (jj = 0; jj < nbin; jj++)
            {
		randx = (xn + (double)rand()/RAND_MAX) * LX;
                randy = (yn + (double)rand()/RAND_MAX) * LY;
                randz = (zn + (double)rand()/RAND_MAX) * LZ;

		for (kk = 0; kk < pointer; kk++)
                {
                   if (fabs(randx - L_atomCoords[0][kk]) <= CUT_OFF)
		      if (fabs(randy - L_atomCoords[1][kk]) <= CUT_OFF)
			 if (fabs(randz - L_atomCoords[2][kk]) <= CUT_OFF) 
                         {
                            sum[tid]++;
                            break;
                         }
                }
 
           }
        }

        for (tid = 0; tid < nthreads; tid++)
        {
//              printf("\ntid: %d, count: %d", tid, sum[tid]);
              totalCount += sum[tid];
        }	

	MPI_Allreduce(&totalCount,&globalCount,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  

        MPI_Finalize();
        t2 = MPI_Wtime();
	if (myid == 0)
	{
        printf("\nPorosity = %lf\n",globalCount/(double)GBIN);
	printf( "Elapsed time is %f\n", t2 - t1 ); 
	}

}
