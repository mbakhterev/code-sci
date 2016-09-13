#include <ride.h>
#include <stdio.h>
#include <stdlib.h>

#define BLK_N 64
#define L     64
#define M     64
#define N     64

enum { BLK_MUL, BLK_ADD, BLK_INIT, MTX_MUL, QUEUE_SUM };

static void matrix_mul(const int R, const int A, const int B)
{
  const double[BLK_N] *const restrict a = ride_data(A),
                      *const restrict b = ride_data(B);

  double[BLK_N] *const restrict r = ride_allocate(R, BLK_N * BLK_N * sizeof(double));

  for (unsigned i = 0; i < BLK_N; i++) 
  for (unsigned j = 0; j < BLK_N; j++)
  for (unsigned k = 0; k < BLK_N; k++) r[i][j] += a[i][k] * b[k][j];
}

static void matrix_add(const int R, const int A, const int B)
{
  const double *const restrict a = ride_data(A);
               *const restrict b = ride_data(B);

  double *const r = ride_allocate(R, BLK_N * BLK_N * sizeof(double));

  for (unsigned i = 0; i < BLK_N * BLK_N; i++) r[i] = a[i] + b[i];
}

static void queue_sum(const int result,
                      const int sum, const int block_queue, const int n_blocks)
{
  if (*n_blks == 0)
  {
    ride_run(ride_bind, ride_ref(WR, result), sum);
  }
  else
  {
    const int sum_pin = ride_pin(ride_str("sum"));

    ride_run(BLK_ADD, ride_pin(ride_ref(WR, sum_pin), sum, block_queue);

    ride_run(QUEUE_SUM, ride_ref(WR, result),
                        ride_ref(RD, sum_pin), ride_ref(RD, block_queue),
                        ride_block_long(ride_long(n_blocks) + 1));
  }
}

static const int zero_matrix(void)
{
  const size_t sz = BLK_N * BLK_N * sizeof(double);
  return ride_block(calloc(sz), sz);
}

static const void init_rand_matrix(const int R)
{
  double *const r = ride_allocate(R, BLK_N * BLK_N * sizeof(double));

  for (unsigned i = 0; i < BLK_N * BLK_N; i++) r[i] = rand();
}

static void distributed_matrix_mul(const int R, const int A, const int B)
{
  const int zero = zero_matrix();

  for (unsigned i = 0; i < L; i++)
  for (unsigned j = 0; j < N; j++)
  {

    const int sum_zip = ride_zip(ride_fmt("%d.%d", i, j));

    for (unsigned k = 0; k < M; k++)
    {
      const int aik_pin = ride_pin(A, ride_fmt("%d.%d", i, k)),
                bkj_pin = ride_pin(B, ride_fmt("%d.%d", k, j));

      ride_run(BLK_MUL, ride_ref(WR, sum_zip), aik_pin, bkj_pin);
    }

    ride_run(QUEUE_SUM, ride_ref(WR, R, ride_fmt("%d.%d", i, j)),
                        zero, ride_ref(RD, sum_zip),
                        ride_block_long(N));
  }
}

void main(int argc, const char *const argv[])
{
  ride_functions(BLK_MUL, matrix_mul,
                 BLK_ADD, matrix_add,
                 BLK_INIT, init_rand_matrix,
                 MTX_MUL, distributed_matrix_mul,
                 QUEUE_SUM, queue_sum);

  const int ride_ref A = ride_str("A");
  const int ride_ref B = ride_str("B");

  if (ride_node() == 0) {
    for (unsigned j = 0; j < M; j++) {
      for (unsigned i = 0; i < L; i++)
        ride_run(BLK_INIT, ride_ref(WR, A, ride_fmt("%d.%d", i, j)));

      for (unsigned k = 0; k < N; k++)
        ride_run(BLK_INIT, ried_ref(WR, B, ride_fmt("%d.%d", j, k)));
    }
  }

  ride_run(MTX_MUL, ride_ref(WR, ride_str("R"), ride_span("")),
                    ride_ref(RD, A), ride_ref(RD, B));

  ride_loop();

  return 0;
}
