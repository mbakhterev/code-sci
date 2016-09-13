#include <ride.h>
#include <stdio.h>

#define BLK_N 64
#define L     64
#define M     64
#define N     64

void matrix_mul(const int R, const int A, const int B)
{
  const double[BLK_N] *const restrict a = ride_data(A),
                      *const restrict b = ride_data(B);

  double[BLK_N] *const restrict r = ride_allocate(R, BLK_N * BLK_N * sizeof(double));

  for (unsigned i = 0; i < BLK_N; i++) 
  for (unsigned j = 0; j < BLK_N; j++)
  for (unsigned k = 0; k < BLK_N; k++) r[i][j] += a[i][k] * b[k][j];
}

void matrix_add(const int R, const int A, const int B)
{
  const double *const restrict a = ride_data(A);
               *const restrict b = ride_data(B);

  double *const r = ride_allocate(R, BLK_N * BLK_N * sizeof(double));

  for (unsigned i = 0; i < BLK_N * BLK_N; i++) r[i] = a[i] + b[i];
}

void queue_sum(const int result,
               const int sum, const int block_queue, const int n_blocks)
{
  if (*n_blks == 0)
  {
    ride_run(ride_bind, ride_ref(WR, result), sum);
  }
  else
  {
    const int sum_pin = ride_pin(ride_str("sum"));

    ride_run(matrix_add, ride_pin(ride_ref(WR, sum_pin), sum, block_queue);

    ride_run(queue_sum, ride_ref(WR, result),
                        ride_ref(RD, sum_pin), ride_ref(RD, block_queue),
                        ride_block_long(ride_long(n_blocks) + 1));
  }
}

const int zero_matrix(void)
{
  const size_t sz = BLK_N * BLK_N * sizeof(double);
  return ride_block(calloc(sz), sz);
}

void distributed_matrix_mul(const int R, const int A, const int B)
{
  const int zero = zero_matrix();

  for (unsigned i = 0; i < L; i++)
  for (unsigned j = 0; i < N; j++)
  {
    const int sum_zip = ride_zip(ride_fmt("%d.%d", i, j));

    for (unsigned k = 0; k < M; k++)
    {
      const int aik_pin = ride_pin(A, ride_fmt("%d.%d", i, k)),
                bkj_pin = ride_pin(B, ride_fmt("%d.%d", k, j));

      ride_run(matrix_mul, ride_ref(WR, sum_zip), aik_pin, bkj_pin);
    }

    ride_run(queue_sum, ride_ref(WR, ride_pin(R, ride_fmt("%d.%d", i, j))),
                        zero, ride_ref(RD, sum_zip),
                        ride_block_long(N));
  }
}
