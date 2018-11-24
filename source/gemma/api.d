// GEMMA C API

module gemma.api;

import bio.std.decompress;
import std.conv;
import std.file;
import std.exception;
import std.stdio;
import std.string;

extern (C++) {

  const __gshared VERSION = "0.10"c;

  import std.algorithm;
  import std.array;
  import std.string : toStringz;

  import std.experimental.logger;
  import core.stdc.string : strncpy;
  import bio.std.range.splitter;

  import gemma.dmatrix;
  import gemma.kinship;
  import gemma.gsl : gsl_matrix;


  // Buf should be preallocated and have at least a size of 16
  char *flmmd_api_version(char *buf) @nogc {
    strncpy(buf,VERSION.ptr,VERSION.length);
    buf[VERSION.length] = '\0';
    return buf;
  }

  void flmmd_compute_bimbam_K(const char *geno_fn, int use_snp_size, int *use_snp, gsl_matrix *k_result) {
    ulong chars = 0;
    ulong token_num = 0;
    double[][] rows;
    rows.length = use_snp_size;
    immutable n_ind = k_result.size1;
    enforce(n_ind == k_result.size2 && n_ind>0,"K matrix result size is invalid");
    info("use snps size ",use_snp_size);
    // ---- Parse the geno file
    info("GZipbyLine");
    auto fn = to!string(fromStringz(cast(char *)geno_fn));
    auto use_snp_num = 0;
    foreach(line, ubyte[] s; GzipbyLine!(ubyte[])(fn)) {
      chars += s.length;
      enforce(use_snp_size >= use_snp_num); // bounds check
      if (use_snp[use_snp_num]) {
        auto tokens = array(SimpleSplitConv!(ubyte[])(s));

        if (token_num == 0) token_num = tokens.length;
        if (token_num != tokens.length) throw new Exception("Number of tokens does not match in line " ~ to!string(line));
        auto elements = new double[token_num-3];
        foreach(i, token; tokens[3..$]) {
          elements[i] = to!double(cast(string)token);
        }
        rows[use_snp_num] = elements;
        use_snp_num++;
      }
    }
    enforce(n_ind == rows[0].length, "Individuals (" ~ to!string(rows[0].length) ~ ") do not match with size of K " ~ to!string(n_ind));
    info("flmmd parsed ",fn," ",use_snp_num," genotypes");
    info("flmmd computes K on ",rows[0].length," individuals");

    // ---- Compute K
    DMatrix G = new DMatrix(rows);
    auto K = kinship_full(G);
  }

  void flmmd_compute_and_write_K(const char *target, const char *geno_fn, const bool is_centered) {
    ulong chars = 0;
    ulong token_num = 0;
    double[][] rows;
    // immutable n_ind = k_result.size1;
    // ---- Parse the geno file
    info("GZipbyLine");
    auto fn = to!string(fromStringz(cast(char *)geno_fn));
    auto use_snp_num = 0;
    foreach(line, ubyte[] s; GzipbyLine!(ubyte[])(fn)) {
      chars += s.length;
      // enforce(use_snp_size >= use_snp_num); // bounds check
      auto tokens = array(SimpleSplitConv!(ubyte[])(s));

      if (token_num == 0) token_num = tokens.length;
      if (token_num != tokens.length) throw new Exception("Number of tokens does not match in line " ~ to!string(line));
      auto elements = new double[token_num-3];
      foreach(i, token; tokens[3..$]) {
        elements[i] = to!double(cast(string)token);
      }
      rows ~= elements;
      use_snp_num++;
    }
    // enforce(n_ind == rows[0].length, "Individuals (" ~ to!string(rows[0].length) ~ ") do not match with size of K " ~ to!string(n_ind));
    info("flmmd parsed ",fn," ",use_snp_num," genotypes");
    info("flmmd computes K on ",rows[0].length," individuals");

    // ---- Compute K
    DMatrix G = new DMatrix(rows);
    auto K = kinship_full(G);

    // ---- Write K

    auto outfn = to!string(fromStringz(cast(char *)target));
    outfn ~= (is_centered ? ".cXX.txt" : ".sXX.txt");
    write_to_file(outfn,cast(immutable DMatrix)K);

  }

} // C++
