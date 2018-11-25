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

  import bio.std.genotype.maf;
  import bio.std.genotype.snp;
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

  /*

    A typical GeneNetwork BIMBAM snp genotype line looks like

  rs31443144, X, Y, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, ...

  */

  void flmmd_compute_and_write_K(const char *target, const char *geno_fn, const char *anno_fn, const bool is_loco, const bool is_centered, const double maf_level) {
    ulong chars = 0;
    ulong token_num = 0;
    double[][] rows;
    // immutable n_ind = k_result.size1;
    // ---- Fetch annotation
    SNP[] snp_annotations;
    if (is_loco)
      snp_annotations = fetch_snp_annotations(to!string(fromStringz(cast(char *)anno_fn)));


    // ---- Parse the geno file
    info("GZipbyLine");
    info("MAF",maf_level);
    auto fn = to!string(fromStringz(cast(char *)geno_fn));
    auto use_snp_num = 0;
    bool[double] genotypes_used;
    foreach(line, ubyte[] s; GzipbyLine!(ubyte[])(fn)) {
      chars += s.length;
      // enforce(use_snp_size >= use_snp_num); // bounds check
      auto tokens = array(SimpleSplitConv!(ubyte[])(s));

      if (token_num == 0) token_num = tokens.length;
      if (token_num != tokens.length) throw new Exception("Number of tokens does not match in line " ~ to!string(line));
      auto elements = new double[token_num-3];
      foreach(i, token; tokens[3..$]) {
        auto t = cast(string)token;
        if (t == "NA") throw new Exception(t); // FIXME - do something
        auto value = to!double(t);
        elements[i] = value;
      }
      // filter on MAF
      auto freqs = maf(elements);
      foreach(k, v; freqs) {
        genotypes_used[k] = true;
      }

      if (maf_level != -1) {
        // Following GEMMA logic
        // FIXME: need to correct for max size
        auto maf = elements.reduce!"a+b"/elements.length;
        if (maf < maf_level || maf > 1.0-maf_level) {
          info(maf,"Dropping ",freqs);
          continue;
        }
        /* but I think this is what it should be
           if (freqs.values.reduce!max > 1.0-default_maf) {
           info("Dropping ",freqs);
           continue;
           }
        */
      }
      // filter on poly
      // FIXME? The code in gemma just checks for multiple types

      // FIXME: CalcHWE

      // FIXME: r2_level

      // Scale matrix
      auto geno_mean = reduce!"a + b"(0.0, elements)/elements.length;
      // writeln(geno_mean);
      auto elements2 = elements.map!(a => a - geno_mean);
      rows ~= array(elements2);
      use_snp_num++;
    }
    info("Genotypes used: ",genotypes_used.keys.sort);
    enforce(is_centered); // FIXME

    // enforce(n_ind == rows[0].length, "Individuals (" ~ to!string(rows[0].length) ~ ") do not match with size of K " ~ to!string(n_ind));
    info("flmmd parsed ",fn," ",use_snp_num," genotypes");
    info("flmmd computes K on ",rows[0].length," individuals");
    if (is_loco) {
      auto chromosomes  = array(snp_annotations.map!(snp => snp.chr)).sort.uniq;
      info("flmmd LOCO on ",chromosomes);
      foreach(chr ; chromosomes) {
        // ---- Compute K
        DMatrix G = new DMatrix(rows);
        auto K = kinship_full(G);

        // ---- Scale K
        foreach (ref e ; K.elements) {
          e *= (1.0/use_snp_num);
        }
      }
    }
    else {
      // ---- Compute K
      DMatrix G = new DMatrix(rows);
      auto K = kinship_full(G);

      // ---- Scale K
      foreach (ref e ; K.elements) {
        e *= (1.0/use_snp_num);
      }
      // ---- Write K

      auto outfn = to!string(fromStringz(cast(char *)target));
      outfn ~= (is_centered ? ".cXX.txt" : ".sXX.txt");
      write_to_file(outfn,cast(immutable DMatrix)K);
    }
  }

} // C++
