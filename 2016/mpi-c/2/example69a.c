// https://parallel.ru/sites/default/files/tech/tech_dev/MPI%26OpenMP/example69a.c

#include <stdio.h>
#include <mpi.h>
#define N 4096
int main(int argc, char* argv[])
{
   int i, j, k, rank, rank_row, rank_col, size, size1, rows, cols, myrows, mycols, last_rows, last_cols, myi, myj, l1, l2, m;
   double t1, t2, sum, s;
   MPI_Comm comm, row, col;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   size1=sqrt(size);
   myrows=mycols=rows=cols=(N-1)/size1+1;
   last_rows=last_cols=N-rows*(size1-1);
   MPI_Comm_split(MPI_COMM_WORLD, rank<size1*size1, rank,
                  &comm);
   if(rank<size1*size1){
      MPI_Comm_split(comm, rank/size1, rank, &row);
      MPI_Comm_split(comm, rank%size1, rank, &col);
      MPI_Comm_rank(row, &rank_row);
      MPI_Comm_rank(col, &rank_col);
      if(rank_row==size1-1) mycols=last_cols;
      if(rank_col==size1-1) myrows=last_rows;
      double **a=(double**)malloc(sizeof(double*)*rows);
      for(i=0; i<rows; i++)
         a[i]=(double*)malloc(sizeof(double)*cols);
      double **b=(double**)malloc(sizeof(double*)*rows);
      for(i=0; i<rows; i++)
         b[i]=(double*)malloc(sizeof(double)*cols);
      double **c=(double**)malloc(sizeof(double*)*rows);
      for(i=0; i<rows; i++)
         c[i]=(double*)malloc(sizeof(double)*cols);
      double *recv_a=(double*)malloc(sizeof(double)*
                                     rows*cols*size1);
      double *recv_b=(double*)malloc(sizeof(double)*
                                     rows*cols*size1);
      for(i=0; i<rows; i++){
         myi=i+rank_col*rows;
         for(j=0; j<cols; j++){
            myj=j+rank_row*cols;
            if(myi<N && myj<N){
               a[i][j]=(myi+myj)*0.0001e0;
               b[i][j]=(myi-myj)*0.0001e0;
            }
            else{
               a[i][j]=0.0e0;
               b[i][j]=0.0e0;
            }
            c[i][j]=0.0e0;
         }
      }
      MPI_Barrier(comm);
      t1=MPI_Wtime();
      l1=rank_row*rows*cols;
      for(i=0; i<rows; i++)
         for(j=0; j<cols; j++)
            recv_a[l1++]=a[i][j];
      l2=rank_col*rows*cols;
      for(i=0; i<rows; i++)
         for(j=0; j<cols; j++)
            recv_b[l2++]=b[i][j];
      for(m=0; m<size1; m++)
         MPI_Bcast(&recv_a[m*rows*cols], rows*cols, MPI_DOUBLE,
                   m, row);
      for(m=0; m<size1; m++)
         MPI_Bcast(&recv_b[m*rows*cols], rows*cols, MPI_DOUBLE,
                   m, col);
      for(m=0; m<size1; m++){
         l1=m*rows*cols;
         for(i=0; i<myrows; i++)
            for(j=0; j<mycols; j++)
               for(k=0; k<cols; k++)
                  c[i][j]+=recv_a[l1+i*cols+k]*
                           recv_b[l1+k*cols+j];
      }
      MPI_Barrier(comm);
      t2=MPI_Wtime();
      for(i=0, sum=0.0e0; i<myrows; i++)
         for(j=0; j<mycols; j++)
            sum+=c[i][j];
      MPI_Reduce(&sum, &s, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
      if(rank==0) printf("N= %d, Nproc=%d, Sum=%lf, Time=%lf\n",
                         N, size, s, t2-t1);
   }
   MPI_Finalize();
}
