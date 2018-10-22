/*
   This code is part of faster_lmm_d and published under the GPLv3
   License (see LICENSE.txt)

   Copyright © 2017 - 2018 Prasun Anand & Pjotr Prins
*/

module gemma.dmatrix;

import std.algorithm;
import std.conv;
import std.math;
import std.stdio;
import std.typecons;

import std.experimental.logger;

alias size_t m_items; // dimensions are never negative

class DMatrix{
  m_items nrows = 0, ncols = 0;
  double[] elements; // 1-dimensional representation
  alias ncols cols;
  alias nrows rows;

  this(m_items num_rows, m_items num_cols) {
    nrows = num_rows;
    ncols = num_cols;
    elements.length = nrows * ncols;
  }

  this(const m_items num_rows, const m_items num_cols, const double[] e) {
    this(num_rows, num_cols);
    elements = cast(double[]) e.dup;
  }

  this(const double[][] rows) {
    this(rows.length,rows[0].length);
    foreach (r, row; rows) {
      foreach (c, col; row) {
        elements[r * ncols + c] = col;
      }
    }
  }

  this(const DMatrix m) {
    this(m.nrows,m.ncols,m.elements);
  }
  pragma(inline) pure const size() { return elements.length; };

}


DMatrix slow_matrix_transpose(const DMatrix m) {
  info("slow_matrix_transpose");
  assert(m.cols > 0 && m.rows > 0);
  if (m.cols == 1 || m.rows == 1) {
    //trace("Fast ",m.cols," x ",m.rows);
    return cast(DMatrix)m;
  }
  else {
    //trace("Slow ", m.cols," x ",m.rows);
    double[] output = new double[m.size];
    auto index = 0;
    auto e = m.elements;
    auto cols = m.cols;
    auto rows = m.rows;
    foreach(i ; 0..cols) {
      foreach(j ; 0..rows) {
        output[index++] = e[j*cols + i];
      }
    }
    return new DMatrix(rows, cols, output);
  }
}
