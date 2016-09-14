#include <ride.h>
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

#define BLK_N 64
#define L     64
#define M     64
#define N     64

void matrix_mul(const int R, const int A, const int B)
{
  const int n = BLK_N;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0,
              ride_data(A), n, ride_data(B), n,
              0.0, ride_allocate(R, n * n * sizeof(double)), n);
}

void matrix_add(const int R, const int A, const int B)
{
  ride_bind_rewrite(R, B);
  cblas_daxpy(BLK_N, 1.0, ride_data(A), 1, ride_data(R), 1);
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

    ride_run("matrix_add", sum_pin, sum, block_queue);

    ride_run("queue_sum", ride_ref(WR, result),
                          ride_ref(RD, sum_pin, block_queue),
                          ride_block_long(ride_long(n_blocks) + 1));
  }
}

static const int zero_matrix(void)
{
  const size_t sz = BLK_N * BLK_N * sizeof(double);
  return ride_block(calloc(sz), sz);
}

void init_rand_matrix(const int R)
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

      ride_run("matrix_mul", ride_ref(WR, sum_zip), aik_pin, bkj_pin);
    }

    ride_run("queue_sum", ride_ref_pin(WR, R, ride_fmt("%d.%d", i, j)),
                          zero, ride_ref(RD, sum_zip),
                          ride_block_long(N));
  }
}

void main(int argc, const char *const argv[])
{
  const int A = ride_str("A");
  const int B = ride_str("B");

  if (ride_node() == 0) {
    for (unsigned j = 0; j < M; j++) {
      for (unsigned i = 0; i < L; i++)
        ride_run("init_rand_matrix", ride_ref_pin(WR, A, ride_fmt("%d.%d", i, j)));

      for (unsigned k = 0; k < N; k++)
        ride_run("init_rand_matrix", ried_ref_pin(WR, B, ride_fmt("%d.%d", j, k)));
    }
  }

  ride_run("distributed_matrix_mul", ride_ref_pin(WR, ride_str("R"), ride_span("")),
                                     ride_ref_pin(RD, ride_seq(A, B));

  ride_loop();

  return 0;
}
