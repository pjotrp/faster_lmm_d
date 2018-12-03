/*
   This code is part of faster_lmm_d and published under the GPLv3
   License (see LICENSE.txt)

   Copyright © 2017 - 2018 Prasun Anand & Pjotr Prins
*/

module gemma.dmatrix;

import std.algorithm;
import std.array;
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
    elements = cast(double[]) e;
  }

  // Copies the contents of rows into elements (FIXME: slow)
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

  ~this() {
    elements = null;
  }

  pragma(inline) pure const size() { return elements.length; };

  auto row(m_items i) const { return elements[i*ncols..(i+1)*ncols]; };

}

struct DMatrixRowIter {
  m_items row = 0;
  immutable DMatrix m;

  this(immutable DMatrix matrix) {
    m = matrix;
  }

  bool empty() const {
    // The range is consumed when begin equals end
    return row == m.rows;
  }

  void popFront() {
    // Skipping the first element is achieved by
    // incrementing the beginning of the range
    ++row;
  }

  const(double[]) front() const {
    // The front element is the one at the beginning
    return m.row(row);
  }

}

DMatrix slow_matrix_transpose(const DMatrix m) {
  info("slow_matrix_transpose");
  assert(m.cols > 0 && m.rows > 0);
  if (m.cols == 1 || m.rows == 1) {
    return cast(DMatrix)m;
  }
  else {
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
    return new DMatrix(cols, rows, output);
  }
}

void write_to_file(string fn, immutable DMatrix m) {
  info("Writing "~fn);
  auto rows = DMatrixRowIter(m);
  auto f = File(fn, "w");
  foreach(ref row; rows) {
    f.writeln(join(row.map!(to!string),"\t"));
  }
}
