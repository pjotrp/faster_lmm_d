/*
   This code is part of faster_lmm_d and published under the GPLv3
   License (see LICENSE.txt)

   Copyright © 2017 - 2018 Prasun Anand & Pjotr Prins
*/

module faster_lmm_d.kinship;

import std.csv;
import std.conv;
import std.exception;
import std.experimental.logger;
import std.file;

import gemma.dmatrix;

DMatrix kinship_full(const DMatrix G)
{
  info("Full kinship matrix used");
  m_items m = G.nrows; // snps
  m_items n = G.ncols; // inds
  log(n," INDs");
  log(m," SNPs");

  DMatrix GT = slow_matrix_transpose(G);
  DMatrix GGT = matrix_mult(GT, G);
  assert(GGT.rows == n);
  assert(GGT.cols == n);

  /*
  info("normalize K");
  DMatrix K = divide_dmatrix_num(mmT, m);

  log("kinship_full K sized ",n," ",K.elements.length);
  log(K.elements[0],",",K.elements[1],",",K.elements[2],"...",K.elements[n-3],",",K.elements[n-2],",",K.elements[n-1]);
  ulong row = n;
  ulong lr = n*n-1;
  ulong ll = (n-1)*n;
  log(K.elements[ll],",",K.elements[ll+1],",",K.elements[ll+2],"...",K.elements[lr-2],",",K.elements[lr-1],",",K.elements[lr]);
  //if(test_kinship){ check_kinship(p_values); }
  check_memory();
  return K;
  */
  return null;
}
