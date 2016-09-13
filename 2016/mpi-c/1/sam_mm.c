// http://www.cs.ucsb.edu/~tyang/class/240b99/homework/sam_mm.c

# include <mpi.h>
# include <stdio.h>
# include <stdlib.h>

static int n;	 	 /* size of the matrix */
static int bsize;	 /* size of the sub-matrix */
static int bnum;	 /* number of sub-matrix columns/rows in the matrix */ 
static int nproc;	 /* number of MPI nodes */
static int myid; 	 /* my own rank */

static double *A, *B, *C, *buf; /* buf is the temporary storage */

int LOCAL(int i)
{
	return i%(n/nproc);
}

int LOCALB(int i)
{
	return i%(bnum/nproc);
}

double *MAP_BLK_A(int i, int j)
{	
	return A+(i*bnum+j)*bsize*bsize;
}

double *MAP_BLK_B(int i, int j)
{	
	return B+(j*bnum+i)*bsize*bsize;
}

double *MAP_BLK_C(int i, int j)
{	
	return C+(i*bnum+j)*bsize*bsize;
}

double *MAP_BLK_BUF(int k)
{	
	return buf+k*bsize*bsize;
}

int OWNER_A(int i, int j)
{
	return i/(n/nproc);
}

int OWNER_B(int i, int j)
{
	return j/(n/nproc);
}

int OWNER_C(int i, int j)
{
	return i/(n/nproc);
}

double * MAP_ELEM_A(int k, int l)
{
	int local_i=LOCAL(k), local_j=l;
	int local_bi=local_i/bsize, local_bj=local_j/bsize;
	int blk_i=local_i%bsize, blk_j=local_j%bsize;
	
	return MAP_BLK_A(local_bi, local_bj)+blk_i*bsize+blk_j;
}

double * MAP_ELEM_B(int k, int l)
{
	int local_i=k, local_j=LOCAL(l);
	int local_bi=local_i/bsize, local_bj=local_j/bsize;
	int blk_i=local_i%bsize, blk_j=local_j%bsize;
	
	return MAP_BLK_B(local_bi, local_bj)+blk_i*bsize+blk_j;
}

double * MAP_ELEM_C(int k, int l)
{
	int local_i=LOCAL(k), local_j=l;
	int local_bi=local_i/bsize, local_bj=local_j/bsize;
	int blk_i=local_i%bsize, blk_j=local_j%bsize;
	
	return MAP_BLK_C(local_bi, local_bj)+blk_i*bsize+blk_j;
}

void do_sub_mm(double *c, double *a, double *b)
{
	int i, j, k;

	for (i=0; i<bsize; i++)
		for (j=0; j<bsize; j++)
			for (k=0; k<bsize; k++)
				c[i*bsize+j]+=a[i*bsize+k]*b[k*bsize+j];
}

void printC()
{
	int i, j;

	for (i=0; i<n; i++)
		if (OWNER_C(i, 0)==myid) {
			for (j=0; j<n; j++)
				printf("%4f  ", *MAP_ELEM_C(i,j));
			printf("\n");
		}
}

int main(int argc, char *argv[])
{
	int i, j, k, token;
	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	if (argc!=3) {
		if (myid==0)
			fprintf(stderr, "Usage, %s <n> <n/q>.\n", argv[0]);
		MPI_Finalize();
		return -1;
	}

	n=atoi(argv[1]);
	bsize=atoi(argv[2]);

	bnum=n/bsize; /* For simplicity, assume they are all correct. */
	if ( (bnum*bsize!=n) || (bnum/nproc*nproc!=bnum) ) {
		if (myid==0)
			fprintf(stderr, "Invalid parameters <nproc=%d> <n=%d> <n/q=%d>.\n", nproc, n, bsize);
		MPI_Finalize();
		return -1;
	}

	A=(double *)malloc(sizeof(double)*(n*n/nproc));
	B=(double *)malloc(sizeof(double)*(n*n/nproc));
	C=(double *)malloc(sizeof(double)*(n*n/nproc));
	buf=(double *)malloc(sizeof(double)*(bsize*n));

	if (!(A && B && C && buf)) {
		if (myid==0)
			fprintf(stderr, "Out of Memory.\n");
		MPI_Finalize();
		return -1;
	}

	memset((char *)C, 0, sizeof(double)*(n*n/nproc));

	if (myid==0) {
		printf("Initialization ...\n");
	}

	for (i=0; i<n; i++)
		for (j=0; j<n; j++) {
			if (OWNER_A(i,j)==myid) {
				*MAP_ELEM_A(i,j)=j;
			}
			if (OWNER_B(i,j)==myid) {
				*MAP_ELEM_B(i,j)=i+j;
			}
	}

	if (myid==0) {
		printf("Done!\n\n");
		printf("Multiplication ...\n");
	}

	for (j=0; j<bnum; j++)
	{
		int root=OWNER_B(0, j*bsize);
		if (root==myid)
			memcpy(buf, MAP_BLK_B(0, LOCALB(j)), sizeof(double)*(bsize*n)) ;

		MPI_Bcast(buf, bsize*n, MPI_DOUBLE, root, MPI_COMM_WORLD);

		for (i=0; i<bnum/nproc; i++) {
			for (k=0; k<bnum; k++) 
				do_sub_mm(MAP_BLK_C(i, j), MAP_BLK_A(i, k), MAP_BLK_BUF(k));
		}
	}

	if (myid==0) {
		printf("Done!\n\n");
	}

	if (n<=10) {
		if (myid==0)
			printf("Printing Matrix C ...\n");

		if (myid!=0) {
			MPI_Recv(&token, 1, MPI_INT, myid-1, 2, MPI_COMM_WORLD, &status);
		}
	
		printC();
		fflush(stdout);
	
		if (myid<nproc-1) {
			MPI_Send(&token, 1, MPI_INT, myid+1, 2, MPI_COMM_WORLD);
		}
	
		if (myid==nproc-1) {
			printf("Done!\n");
		}
	}

	MPI_Finalize();

	return 0;
}
