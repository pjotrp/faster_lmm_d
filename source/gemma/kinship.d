/*
   This code is part of faster_lmm_d and published under the GPLv3
   License (see LICENSE.txt)

   Copyright © 2017 - 2018 Prasun Anand & Pjotr Prins
*/

module gemma.kinship;

import std.array : array;
import std.algorithm;
import std.csv;
import std.conv;
import std.exception;
import std.experimental.logger;
import std.file;

import std.experimental.logger;

import gemma.blas;
import gemma.dmatrix;

DMatrix compute_K(const DMatrix G)
{
  info("Full kinship matrix used");
  m_items m = G.nrows; // snps
  m_items n = G.ncols; // inds
  log(n," INDs");
  log(m," SNPs");

  DMatrix GT = slow_matrix_transpose(G);
  assert(GT.rows == G.cols);
  assert(GT.cols == G.rows);
  DMatrix K = matrix_mult(GT, G);
  info("DONE rows is ",K.rows," cols ",K.cols);
  assert(K.rows == n);
  assert(K.cols == n);

  // ---- Scale K
  foreach (ref e ; K.elements) {
    e *= (1.0/G.rows);
  }

  return K;
}



import gemma.api : SnpGenotypes;

DMatrix compute_K(const SnpGenotypes[] rows, const string chr) {
  DMatrix G = new DMatrix(array(rows.filter!(r => r.snp.chr != chr ).map!(r => r.genotypes)));
  return compute_K(G);
}
