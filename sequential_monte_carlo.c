//#pragma warning(disable:4996)
#include "mpi.h"
#include<stdio.h>
#include <stdlib.h>
#include<math.h>

#define CUT_OFF 3.25
#define X 64
#define Y 64
#define Z 64
#define LEN(x)  (sizeof(x) / sizeof((x)[0]))

struct AtomData {
	int noOfAtoms;
	float* atomCoords[4];
}AtomDataObj;
//typedef struct AtomData AtomData;

/*long getHits(int noOfRandomPts) {
	int hitsCount = 0;
	float x, y, z;
	FILE *fop = fopen("mc.out", "w");
	srand((unsigned int)time(NULL));
	for (int i = 0; i < noOfRandomPts; i++) {
		//generate point
		x = ((double)rand() / (double)(RAND_MAX)) * X;
		y = ((double)rand() / (double)(RAND_MAX)) * Y;
		z = ((double)rand() / (double)(RAND_MAX)) * Z;

		//check if within sys or outside
		for (int j = 0; j < AtomDataObj.noOfAtoms; j++) {
			float xa = AtomDataObj.atomCoords[0][j];
			float ya = AtomDataObj.atomCoords[1][j];
			float za = AtomDataObj.atomCoords[2][j];
			if ((abs(x - xa) <= CUT_OFF) && (abs(y - ya) <= CUT_OFF) && (abs(z - za) <= CUT_OFF)) {
				printf("\natom coords %f %f %f", xa, ya, za);
				printf("\nrandom Pt %f %f %f", x, y, z);
				hitsCount++;
				break;
			}
		}

	}
	return hitsCount;

}*/

/*long getHitsBySequentialPoints() {
	int hitsCount = 0;
	float xyz[3];
	FILE *fop = fopen("mc.out", "w");
	//to start with different seed for random generator
	//srand((unsigned int)time(NULL));
	//generate point
	int x,y,z,i;
	for (x = 0; x < X;) {
		for ( y = 0; y < Y;) {
			for ( z = 0; z < Z;) {
				//check if within sys or outside
				for ( i = 0; i < AtomDataObj.noOfAtoms; i++) {
					float xa = AtomDataObj.atomCoords[0][i];
					float ya = AtomDataObj.atomCoords[1][i];
					float za = AtomDataObj.atomCoords[2][i];
					if ((abs(x - xa) <= CUT_OFF) && (abs(y - ya) <= CUT_OFF) && (abs(z - za) <= CUT_OFF)) {
						//fprintf(fop, "\natom coords %f %f %f", xa, ya, za);
						//fprintf(fop, "\nrandom Pt %f %f %f", xyz[0], xyz[1], xyz[2]);
						printf("\natom coords %f %f %f", xa, ya, za);

						hitsCount++;
						break;
					}
				}
				z += CUT_OFF;

			}
			y += CUT_OFF;
		}
		x += CUT_OFF;

	}
	return hitsCount;

}*/
long getHitsByPSequentialPoints(int procNo) {
	int hitsCount = 0;
	float xyz[3];
	FILE *fop = fopen("mc.out", "w");
	//
	int procs[3] = { 2,2,2 };
	int startX = procNo / (procs[1] * procs[2]);
	startX *= X / procs[0];

	int startY = (procNo / procs[2]) % procs[1];
	startY *= Y / procs[1];
	int startZ = procNo % procs[2];
	startZ *= Z / procs[2];
	//	printf("\n%d,%d,%d", startX, startY, startZ);

		//generate point
	float x, y, z;
	int i;
	for (x = startX; x < startX + X / procs[0];) {
		for (y = startY; y < startY + Y / procs[1];) {
			for (z = startZ; z < startZ + Z / procs[2];) {
				//check if within sys or outside
				for (i = 0; i < AtomDataObj.noOfAtoms; i++) {
					float xa = AtomDataObj.atomCoords[0][i];
					float ya = AtomDataObj.atomCoords[1][i];
					float za = AtomDataObj.atomCoords[2][i];
					if ((abs(x - xa) <= CUT_OFF) && (abs(y - ya) <= CUT_OFF) && (abs(z - za) <= CUT_OFF)) {
						//fprintf(fop, "\natom coords %f %f %f", xa, ya, za);
						//fprintf(fop, "\nrandom Pt %f %f %f", xyz[0], xyz[1], xyz[2]);
						//printf("\natom coords %f %f %f", xa, ya, za);

						hitsCount++;
						break;
					}
				}
				z += CUT_OFF;

			}
			y += CUT_OFF;
		}
		x += CUT_OFF;

	}
	return hitsCount;

}


void readAtoms() {
	/* Open an MD-configuration file */
	FILE *fp = fopen("atom_coords.in", "r");
	int i;
	/* Read the # of atoms */
	fscanf(fp, "%d", &AtomDataObj.noOfAtoms);
	for (i = 0; i < 4; i++) {
		AtomDataObj.atomCoords[i] = (float*)malloc(AtomDataObj.noOfAtoms * sizeof(float));
	}
	for (i = 0; i < AtomDataObj.noOfAtoms; i++) {
		fscanf(fp, "%f%f%f%f", &AtomDataObj.atomCoords[0][i], &AtomDataObj.atomCoords[1][i], &AtomDataObj.atomCoords[2][i], &AtomDataObj.atomCoords[3][i]);

	}
	//printf("%f", AtomDataObj.atomCoords[0][2]);
}

double global_sum(double partial, int myid, int nprocs) {
	/* Implement your own global summation here */
	double mydone = partial, hisdone;
	int bitValue;
	MPI_Status status;
	for (bitValue = 1; bitValue < nprocs; bitValue *= 2) {
		int partner = myid ^ bitValue;
		MPI_Send(&mydone, 1, MPI_DOUBLE, partner, bitValue, MPI_COMM_WORLD);
		MPI_Recv(&hisdone, 1, MPI_DOUBLE, partner, bitValue, MPI_COMM_WORLD, &status);
		mydone += hisdone;
	}
	return mydone;
}

int main(int argc, char *argv[]) {
	readAtoms();
	//printf("Hello Vis");
	int noOfRandomPts = 500000;
	int nprocs;  /* Number of processors */
	int myid;    /* My rank */

	//long hits=getHits(noOfRandomPts);
	//printf("\n%f", (float)(1.0* hits/ noOfRandomPts));
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	long hits = getHitsByPSequentialPoints(myid);
	int totalHits = global_sum(hits, myid, nprocs);
	printf("\nhits-%d", hits);

	if (myid == 0) {
		long totalPoints = (X / CUT_OFF)*(Y / CUT_OFF)*(Z / CUT_OFF);
		printf("\ntotalHits = %d", totalHits);
		printf("\nTotalPoints = %f", (float)(1.0* hits / (totalPoints)));
	}

	MPI_Finalize();

	return 0;

	/*long totalPoints = (X / CUT_OFF)*(Y / CUT_OFF)*(Z / CUT_OFF);
	printf("\n%f", (float)(1.0* hits / (totalPoints)));*/
	//printf("\n%f", (float)(1.0* hits/ noOfRandomPts));


}
