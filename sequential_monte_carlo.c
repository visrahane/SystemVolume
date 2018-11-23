//#pragma warning(disable:4996)
#include "mpi.h"
#include<stdio.h>
#include <stdlib.h>
#include <math.h>

#define CUT_OFF 1.5
#define X 1060.28
#define Y 212.06
#define Z 212.06
#define LEN(x)  (sizeof(x) / sizeof((x)[0]))
long totalPoints = 0;

struct AtomData {
	int noOfAtoms;
	float* atomCoords[3];
}AtomDataObj,ProcAtomObj;

long getHitsByPSequentialPoints(int procNo) {
	int hitsCount = 0;
	float xyz[3];
	//	FILE *fop = fopen("mc.out", "w");
	//
	int procs[3] = { 20,2,2 };
	int startX = procNo / (procs[1] * procs[2]);
	float XPerProc = X / procs[0];
	startX *= XPerProc;

	int startY = (procNo / procs[2]) % procs[1];
	float YPerProc = Y / procs[1];
	startY *= YPerProc;

	int startZ = procNo % procs[2];
	float ZPerProc = Z / procs[2];
	startZ *= ZPerProc;

	int countOfAtoms=fetchLocalAtoms(procs, procNo);
	//generate point
	float x, y, z;
	int i, points = 0;
	double t1 = MPI_Wtime();
		
	for (x = startX; x < startX + XPerProc; x += CUT_OFF) {
		for (y = startY; y < startY + YPerProc; y += CUT_OFF) {
			for (z = startZ; z < startZ + ZPerProc; z += CUT_OFF) {
				//check if within sys or outside
				points++;
				for (i = 0; i < countOfAtoms; i++) {
					float xa = ProcAtomObj.atomCoords[0][i];
					float ya = ProcAtomObj.atomCoords[1][i];
					float za = ProcAtomObj.atomCoords[2][i];
					if ((fabs(x - xa) <= CUT_OFF) && (fabs(y - ya) <= CUT_OFF) && (fabs(z - za) <= CUT_OFF)) {
						//fprintf(fop, "\natom coords %f %f %f", xa, ya, za);
						//fprintf(fop, "\nrandom Pt %f %f %f", xyz[0], xyz[1], xyz[2]);
						//printf("\natom coords %f %f %f", xa, ya, za);

						hitsCount++;
						break;
					}
				}
			}
		}
	}
	printf("time for for loop is %f\n", MPI_Wtime() - t1);
	printf("points:%d\n", points);
	MPI_Allreduce(&points, &totalPoints, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	return hitsCount;

}

int fetchLocalAtoms(int procs[],int myid) {
	int i, count=0,xn,yn,zn,sid;
	float XPerProc = X / procs[0];
	float YPerProc = Y / procs[1];
	float ZPerProc = Z / procs[2];
	for (i = 0; i < 3; i++)
	{
		ProcAtomObj.atomCoords[i] = (float*)malloc(AtomDataObj.noOfAtoms * sizeof(float));
	}

	for (i = 0; i < AtomDataObj.noOfAtoms; i++)
	{
		xn = (int)(AtomDataObj.atomCoords[0][i] / XPerProc);
		yn = (int)(AtomDataObj.atomCoords[1][i] / YPerProc);
		zn = (int)(AtomDataObj.atomCoords[2][i] / ZPerProc);
		sid = xn * procs[1] * procs[2] + yn * procs[2] + zn;
		if (sid == myid)
		{
			ProcAtomObj.atomCoords[0][count] = AtomDataObj.atomCoords[0][i];
			ProcAtomObj.atomCoords[1][count] = AtomDataObj.atomCoords[1][i];
			ProcAtomObj.atomCoords[2][count] = AtomDataObj.atomCoords[2][i];
			count++;
		}
	}
	return count;
}
void readAtoms(char *file_name) {
	/* Open an MD-configuration file */
	FILE *fp = fopen(file_name, "r");
	int i;
	/* Read the # of atoms */
	fscanf(fp, "%d", &AtomDataObj.noOfAtoms);
	for (i = 0; i < 3; i++) {
		AtomDataObj.atomCoords[i] = (float*)malloc(AtomDataObj.noOfAtoms * sizeof(float));
	}
	for (i = 0; i < AtomDataObj.noOfAtoms; i++) {
		fscanf(fp, "%f%f%f", &AtomDataObj.atomCoords[0][i], &AtomDataObj.atomCoords[1][i], &AtomDataObj.atomCoords[2][i]);

	}

}

int main(int argc, char *argv[]) {
	double t1 = MPI_Wtime();
	readAtoms(argv[1]);

	int noOfRandomPts = 500000;
	int nprocs;  /* Number of processors */
	int myid;    /* My rank */

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	long hits = getHitsByPSequentialPoints(myid);
	int totalHits;
	MPI_Allreduce(&hits, &totalHits, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	double t2 = MPI_Wtime();
	if (myid == 0) {
		//long totalPoints = (X / (2 * CUT_OFF))*(Y / (2 * CUT_OFF))*(Z / (2 * CUT_OFF));
		printf("\ntotalHits = %d\n", totalHits);
		printf("\ntotalPoints = %d\n", totalPoints);
		printf("\nVolume = %f\n", (float)(1.0* totalHits / (totalPoints)));
		printf("Total Elapsed time is %f\n", t2 - t1);

	}
	MPI_Finalize();

	return 0;
}

