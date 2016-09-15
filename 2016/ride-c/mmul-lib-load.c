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
              rval(A), n, rval(B), n, 0.0, ralloc(R, n*n * sizeof(double)), n);
}

void matrix_add(const int R, const int A, const int B)
{
  rbind(R, B);
  cblas_daxpy(BLK_N, 1.0, rval(A), 1, rval(R), 1);
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

static void distributed_matrix_mul(const int R, const int A, const int B)
{
  const int Z = zero_matrix();

  for (unsigned i = 0; i < L; i++)
  for (unsigned j = 0; j < N; j++) {
    const int sum_zip = rzip(rfmt("%d.%d", i, j));

    for (unsigned k = 0; k < M; k++)
    {
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
  const int A = rstr(argc > 2 ? argv[2] : "A"),
            B = rstr(argc > 3 ? argv[3] : "B"),
            R = rstr(argc > 1 ? argv[1] : "R");

  if (ride_node() == 0) {
    ride_run("distributed_matrix_mul", rrefpin(WR, R, ride_span("")),
                                       rrefpin(RD, RROOT, rseq(A, B)));

    rrun("bind", rref(WR, rzip(RAPI, rstr("update-root"))), rpinref(RD));
  }

  return 0;
}
