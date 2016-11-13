#include <ride.h>
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

static const int BLK_N = 64, L = 64, M = 64, N = 64;

void matrix_mul(const int R, const int A, const int B)
{
  const int n = BLK_N;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0,
              rval(A), n, rval(B), n, 0.0, ralloc(R, n*n * sizeof(double)), n);
}

void matrix_add(const int R, const int A, const int B)
{
  rbind(R, A);
  cblas_daxpy(BLK_N, 1.0, rval(B), 1, rval(R), 1);
}

void queue_sum(const int result,
               const int sum, const int block_queue, const int n_blocks)
{
  if (n_blocks == 0) rrun("bind", rref(WR, result), sum);
  else {
    const int sum_pin = rpin(rstr("sum"));

    rrun("matrix_add", sum_pin, sum, block_queue);

    rrun("queue_sum", rref(WR, result), rref(RD, sum_pin, block_queue),
                      rblocklong(rvallong(n_blocks) + 1));
  }
}

static const int zero_matrix(void)
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

      ride_run("matrix_mul", rref(WR, sum_zip), aik_pin, bkj_pin);
    }

    rrun("queue_sum", rrefpin(WR, R, ride_fmt("%d.%d", i, j)),
                              Z, rref(RD, sum_zip), rblocklong(N));
  }
}

void watchdog(const int proc_zip, const int start)
{
  if(rvallong(proc_zip) != RDONE) rrun("watchdog", proc_zip, start);
  else
    printf("time: %d\n", difftime(time(NULL), rvallong(start)));
}

void main(int argc, const char *const argv[])
{
  rrun("watchdog", rzip(rstr("proc")), rblocklong(time(NULL)));

  for (unsigned j = 0; j < M; j++) {
    for (unsigned i = 0; i < L; i++)
      rrun("init_rand_matrix", rrefpin(WR, rstr("A"), rfmt("%d.%d", i, j)));

    for (unsigned k = 0; k < N; k++)
      rrun("init_rand_matrix", rrefpin(WR, rstr("B"), rfmt("%d.%d", j, k)));
  }

  rrun("fork", rrefzip(WR, rstr("proc")),
       "distributed_matrix_mul", rrefpin(WR, rstr("R"), rspan("")),
                                 rrefpin(RD, rbr(rstr("A"), rstr("B"))));
  return 0;
}
