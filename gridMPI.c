#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>



void main(int argc,char **argv) 
{

	double t1, t2;
	int nAtoms, myid, nproc, totalCount = 0, globalCount, accum = 0, globalAccum;
	int i;
	int NX = 2, NY = 2, NZ = 2;
//	double X = 1060.28, Y = 212.06, Z = 212.06;
//      double X = 1053.02, Y = 210.62, Z = 210.62;
	double X = 54.26, Y = 54.26, Z =54.26;
	double CUT_OFF = 1.38;
	double STEP = 1.38;
	double LX = X/NX;
	double LY = Y/NY;
        double LZ = Z/NZ;	
        double* atomCoords[3];


        FILE *fp = fopen("cuzr.txt", "r");
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
        t1 = MPI_Wtime();

	int ii, kk;
        int xn, yn, zn, sid, pointer=0;
        double ixx, iyy, izz;
	double* L_atomCoords[3];

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
		if (sid == myid)
		 {
		    L_atomCoords[0][pointer] = atomCoords[0][ii];
                    L_atomCoords[1][pointer] = atomCoords[1][ii];               
		    L_atomCoords[2][pointer] = atomCoords[2][ii];
		    pointer ++;
		 }
	    }		
	//    printf("POINTER: %d\n", pointer);
            xn = myid / (NY*NZ);
            yn = (myid / NY) % NZ;
            zn = myid % NZ;
	    for (ixx = xn*LX; ixx < xn*LX + LX; ixx += STEP)
            for (iyy = yn*LY; iyy < yn*LY + LY; iyy += STEP)
            for (izz = zn*LZ; izz < zn*LZ + LZ; izz += STEP)
		{
		   accum++;
	    	   for (kk = 0; kk < pointer; kk++)
                   {
		      if (fabs(ixx - L_atomCoords[0][kk]) <= CUT_OFF)
		      if (fabs(iyy - L_atomCoords[1][kk]) <= CUT_OFF)
	              if (fabs(izz - L_atomCoords[2][kk]) <= CUT_OFF) 
                         {
                            totalCount++;
                            break;
                         }
                   }
		}

	MPI_Allreduce(&totalCount,&globalCount,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  
        MPI_Allreduce(&accum,&globalAccum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        t2 = MPI_Wtime();

	if (myid == 1)
	{
           printf("\nPorosity = %lf\n", globalCount/(double)globalAccum);
	   printf( "Elapsed time is %f\n", t2 - t1 ); 
	}

        MPI_Finalize();
}
