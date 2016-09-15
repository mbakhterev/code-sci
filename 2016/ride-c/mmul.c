#include <ride.h>
#include <stdio.h>
#include <stdlib.h>

#define BLK_N 64
#define L     64
#define M     64
#define N     64

void matrix_mul(const int R, const int A, const int B)
{
  const double (*const restrict a)[BLK_N] = rval(A),
               (*const restrict b)[BLK_N] = rval(B);

  double (*const restrict r)[BLK_N] = ralloc(R, BLK_N * BLK_N * sizeof(double));

  for (unsigned i = 0; i < BLK_N; i++) 
  for (unsigned j = 0; j < BLK_N; j++)
  for (unsigned k = 0; k < BLK_N; k++) r[i][j] += a[i][k] * b[k][j];
}

void matrix_add(const int R, const int A, const int B)
{
  const double *const restrict a = rval(A);
               *const restrict b = rval(B);

  double *const r = ralloc(R, BLK_N * BLK_N * sizeof(double));

  for (unsigned i = 0; i < BLK_N * BLK_N; i++) r[i] = a[i] + b[i];
}

void queue_sum(const int result,
               const int sum, const int block_queue, const int n_blocks)
{
  if (n_blocks == 0) rrun("bind", rref(WR, result), sum);
  else {
    const int sum_pin = rpin(ride_str("sum"));

    rrun("matrix_add", rref(sum_pin), sum, block_queue);

    rrun("queue_sum", rref(WR, result), rref(RD, sum_pin, block_queue),
                      rblocklong(rvallong(n_blocks) + 1));
  }
}

static int zero_matrix(void)
{
  const size_t sz = BLK_N * BLK_N * sizeof(double);
  return rblock(calloc(sz), sz);
}

void init_rand_matrix(const int R)
{
  double *const r = ralloc(R, BLK_N * BLK_N * sizeof(double));

  for (unsigned i = 0; i < BLK_N * BLK_N; i++) r[i] = rand();
}

static void distributed_matrix_mul(const int R, const int A, const int B)
{
  const int Z = zero_matrix();

  for (unsigned i = 0; i < L; i++)
  for (unsigned j = 0; j < N; j++) {
    const int sum_zip = rzip(rfmt("%d.%d", i, j));

    for (unsigned k = 0; k < M; k++) {
      const int aik_pin = rpin(A, rfmt("%d.%d", i, k)),
                bkj_pin = rpin(B, rfmt("%d.%d", k, j));

      rrun("matrix_mul", rref(WR, sum_zip), aik_pin, bkj_pin);
    }

    rrun("queue_sum", rrefpin(WR, R, rfmt("%d.%d", i, j)),
                      Z, rref(RD, sum_zip), rblocklong(N));
  }
}

void main(int argc, const char *const argv[])
{
  if (rnode() == 0) {
    for (unsigned j = 0; j < M; j++) {
      for (unsigned i = 0; i < L; i++)
        rrun("init_rand_matrix", rpin(WR, rstr("A"), rfmt("%d.%d", i, j)));

      for (unsigned k = 0; k < N; k++)
        rrun("init_rand_matrix", rpin(WR, rstr("B"), rfmt("%d.%d", j, k)));
    }

    rrun("distributed_matrix_mul", rrefpin(WR, rstr("R"), rspan("")),
                                   rrefpin(RD, rseq(rstr("A"), rstr("B"))));
  }

  return 0;
}
