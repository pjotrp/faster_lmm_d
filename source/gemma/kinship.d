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
import std.stdio;

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

  DMatrix K = matrix_mult_transpose(G);
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
  auto select = new bool[rows.length];
  if (chr == "all") {
    select[] = true;
  }
  else
    foreach (i, row ; rows) {
      select[i] = (row.snp.chr != chr);
    }
  auto count = reduce!"a + to!int(b)"(0,select);
  info("Counted ",count," genotypes for ",chr);
  DMatrix G = new DMatrix(count,rows[0].genotypes.length);
  auto r = 0;
  foreach(i, selected ; select) {
    if (selected) {
      auto row = rows[i];
      foreach(c, geno; row.genotypes) {
        G.elements[r * G.ncols + c] = geno;
      }
      r++;
    }
  }
  auto K = compute_K(G);
  G.elements = null; // tell GC to free up
  return K;
}
