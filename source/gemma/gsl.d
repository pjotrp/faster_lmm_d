/*
   This code is part of faster_lmm_d and published under the GPLv3
   License (see LICENSE.txt)

   Copyright © 2017 - 2018 Prasun Anand & Pjotr Prins
*/

import std.conv;
import std.experimental.logger;

extern (C) {

  struct gsl_block
  {
    size_t size;
    double *data;
  };

  // see https://www.gnu.org/software/gsl/manual/html_node/Matrices.html#Matrices
  struct gsl_matrix
  {
    size_t size1; // rows
    size_t size2; // cols
    size_t tda;   // trailing dimension (actual row length including trailing places)
    double *data;
    gsl_block *block;
    int owner;    // block is owned by gsl_matrix
  };

}
