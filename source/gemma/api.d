// GEMMA C API

module gemma.api;

import bio.std.decompress;
import std.conv;
import std.exception;
import std.stdio;

extern (C) {

  import std.array;
  import std.string : toStringz;
  import core.stdc.string : strncpy;
  import bio.std.range.splitter;

  import gemma.dmatrix;
  import gemma.kinship;

  const VERSION = "0.10"c;

  // Buf should be preallocated and have at least a size of 16
  char *flmmd_api_version(char *buf) @nogc {
    strncpy(buf,VERSION.ptr,VERSION.length);
    buf[VERSION.length] = '\0';
    return buf;
  }

  void flmmd_compute_bimbam_K() {
    ulong lines = 0;
    ulong chars = 0;
    ulong token_num = 0;
    double[][] rows;
    foreach(ubyte[] s; GzipbyLine!(ubyte[])("example/mouse_hs1940.geno.txt.gz")) {
      // test file contains 7320 lines 4707218 characters
      // write(cast(string)s);
      chars += s.length;
      lines += 1;
      // writeln(s);
      auto tokens = array(SimpleSplitConv!(ubyte[])(s));

      if (token_num == 0) token_num = tokens.length;
      if (token_num != tokens.length) throw new Exception("Number on tokens does not match in line " ~ to!string(lines));

      auto elements = new double[token_num-3];
      foreach(i, token; tokens[3..$]) {
        elements[i] = to!double(cast(string)token);
      }
      rows ~= elements;
    }
    DMatrix G = new DMatrix(rows);
    auto K = kinship_full(G);

    assert(chars == 71434128,"chars " ~ to!string(chars));
    assert(lines == 12226,"lines " ~ to!string(lines));
    assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"hello, 1 2 \n\t3  4 \n")) == ["hello","1","2","3","4"]);
  }


} // extern C
