// GEMMA C API

module gemma.api;

import bio.std.decompress;
import std.conv;
import std.exception;
import std.stdio;
import std.string;

extern (C) {
  const VERSION = "0.10"c;
}

extern (C++) {

  import std.algorithm;
  import std.array;
  import std.string : toStringz;

  import std.experimental.logger;
  import core.stdc.string : strncpy;
  import bio.std.range.splitter;

  import gemma.dmatrix;
  import gemma.kinship;


  // Buf should be preallocated and have at least a size of 16
  char *flmmd_api_version(char *buf) @nogc {
    strncpy(buf,VERSION.ptr,VERSION.length);
    buf[VERSION.length] = '\0';
    return buf;
  }

  void flmmd_compute_bimbam_K(const char *geno_fn, int use_snp_size, int *use_snp) {
    ulong lines = 0;
    ulong chars = 0;
    ulong token_num = 0;
    double[][] rows;
    info("use snps size ",use_snp_size);
    info("GZipbyLine");
    auto fn = to!string(fromStringz(cast(char *)geno_fn));
    foreach(ubyte[] s; GzipbyLine!(ubyte[])(fn)) {
      // write(cast(string)s);

      chars += s.length;
      lines += 1;
      // writeln(s);
      enforce(use_snp_size >= lines);
      if (use_snp[lines]) {
        auto tokens = array(SimpleSplitConv!(ubyte[])(s));

        if (token_num == 0) token_num = tokens.length;
        if (token_num != tokens.length) throw new Exception("Number on tokens does not match in line " ~ to!string(lines));

        auto elements = new double[token_num-3];
        foreach(i, token; tokens[3..$]) {
          elements[i] = to!double(cast(string)token);
        }
        rows ~= elements;
      }
      // writeln("");
    }
    DMatrix G = new DMatrix(rows);

    info("flmmd parsed ",fn," ",lines," genotypes");
    info("flmmd computes K on ",rows[0].length," individuals");
    auto K = kinship_full(G);

    // assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"hello, 1 2 \n\t3  4 \n")) == ["hello","1","2","3","4"]);

    // writeln(G.row(0));
    // writeln(reduce!("a + b")(G.row(0)));
    // assert(to!int(reduce!("a + b")(G.row(0))) == 1721 );
  }


} // extern C
