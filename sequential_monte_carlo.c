#pragma warning(disable:4996)
#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#define L 3.25;
#define X 64;
#define Y 64;
#define Z 64;
#define LEN(x)  (sizeof(x) / sizeof((x)[0]))

struct AtomData {
	int noOfAtoms;
	float* atomCoords[4];
}AtomDataObj;
//typedef struct AtomData AtomData;

long getHits(int noOfRandomPts) {
	int hitsCount = 0;
	float x, y, z;
	FILE *fop = fopen("mc.out", "w");
	srand((unsigned int)time(NULL));
	for (int i = 0; i < noOfRandomPts; i++) {
		//generate point
		x = ((double)rand() / (double)(RAND_MAX)) * X ;
		y = ((double)rand() / (double)(RAND_MAX)) * Y;
		z = ((double)rand() / (double)(RAND_MAX)) * Z;
		
		//check if within sys or outside
		for (int i = 0; i < AtomDataObj.noOfAtoms; i++) {
			float xa = AtomDataObj.atomCoords[0][i];
			float ya = AtomDataObj.atomCoords[1][i];
			float za = AtomDataObj.atomCoords[2][i];
			if ((abs(x - xa) <= 3.25) && (abs(y - ya) <= 3.25) && (abs(z - za) <= 3.25)) {
				printf("\natom coords %f %f %f", xa,ya,za);
				printf("\nrandom Pt %f %f %f", x, y, z);
				hitsCount++;
				break;
			}
		}

	}
	return hitsCount;

}
void readAtoms() {
	/* Open an MD-configuration file */
	FILE *fp = fopen("atom_coords.in", "r");
	
	/* Read the # of atoms */
	fscanf(fp, "%d", &AtomDataObj.noOfAtoms);
	for (int i = 0; i < 4; i++) {
		AtomDataObj.atomCoords[i] = (float*)malloc(AtomDataObj.noOfAtoms*sizeof(float));
	}
	for (int i = 0; i < AtomDataObj.noOfAtoms; i++) {
		fscanf(fp, "%f%f%f%f", &AtomDataObj.atomCoords[0][i], &AtomDataObj.atomCoords[1][i], &AtomDataObj.atomCoords[2][i], &AtomDataObj.atomCoords[3][i]);

	}
	//printf("%f", AtomDataObj.atomCoords[0][2]);
}

void main() {
	readAtoms();
	//printf("Hello Vis");
	int noOfRandomPts = 500000;
	long hits=getHits(noOfRandomPts);
	printf("\n%f", (float)(1.0* hits/ noOfRandomPts));

}