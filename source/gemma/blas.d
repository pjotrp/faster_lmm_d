/*
   This code is part of faster_lmm_d and published under the GPLv3
   License (see LICENSE.txt)

   Copyright © 2017 - 2018 Prasun Anand & Pjotr Prins
*/

import std.conv;

import std.experimental.logger;

import gemma.dmatrix;
// import cblas : gemm;

enum CBLAS_ORDER : blasint  {
  ///
  RowMajor=101,
  ///
  ColMajor=102
};

///
alias Order = CBLAS_ORDER;

///
enum CBLAS_TRANSPOSE : blasint {
  ///
  NoTrans=111,
  ///
  Trans=112,
  ///
  ConjTrans=113,
  ///
  ConjNoTrans=114
};

///
alias Transpose = CBLAS_TRANSPOSE;

alias blasint = int;

extern (C) {
void cblas_dgemm(in CBLAS_ORDER order, in CBLAS_TRANSPOSE TransA, in CBLAS_TRANSPOSE TransB, in blasint M, in blasint N, in blasint K, in double alpha, in double *A, in blasint lda, in double *B, in blasint ldb, in double beta, double *C, in blasint ldc);
}

DMatrix matrix_mult(const DMatrix lha,const DMatrix rha) {
  info("matrix_mult");
  double[] C = new double[lha.nrows*rha.ncols];
  cblas_dgemm(Order.RowMajor, Transpose.NoTrans, Transpose.NoTrans, to!int(lha.rows), to!int(rha.cols), to!int(lha.cols), /*no scaling*/
              1,lha.elements.ptr, to!int(lha.ncols), rha.elements.ptr, to!int(rha.ncols), /*no addition*/0, C.ptr, to!int(rha.ncols));
  // auto res_shape = [lha.rows(),rha.cols()];
  return new DMatrix(lha.rows, rha.cols, C);
}


/*
void matrixtTransposeMult(std::vector<float>& inMatrix, int rowNum)
{
    cublasHandle_t handle;
    cublasCreate(&handle);

    int colNum = (int)inMatrix.size() / rowNum;
    thrust::device_vector<float> d_InMatrix(inMatrix);
    thrust::device_vector<float> d_outputMatrix(rowNum*rowNum);
    float alpha = 1.0f;
    float beta = 0.0f;

    cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, rowNum, rowNum, colNum, &alpha,
        thrust::raw_pointer_cast(d_InMatrix.data()), rowNum, thrust::raw_pointer_cast(d_InMatrix.data()), rowNum, &beta,
        thrust::raw_pointer_cast(d_outputMatrix.data()), rowNum);

    thrust::host_vector<float> result = d_outputMatrix;
    for (auto elem : result)
        std::cout << elem << ",";
    std::cout << std::endl;

    cublasDestroy(handle);
}
*/

DMatrix matrix_mult_transpose(const DMatrix G) {
  info("matrix_mult_transpose");
  int rows = to!int(G.cols);
  int cols = to!int(G.rows);
  double[] C = new double[rows*rows]; // result matrix
  cblas_dgemm(Order.RowMajor, Transpose.Trans, Transpose.NoTrans, rows, rows, cols, /*no scaling*/
              1.0, G.elements.ptr, rows, G.elements.ptr, rows, /*no addition*/0.0, C.ptr, rows);
  return new DMatrix(rows,rows,C);
}
