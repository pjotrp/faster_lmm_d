// GEMMA C API

module gemma.api;

import bio.std.decompress;

// import std.concurrency;
import std.conv;
import std.file;
import std.exception;
import std.parallelism;
import std.stdio;
import std.string;
import std.typecons;

alias SnpGenotypes = Tuple!(SNP, "snp", double[], "genotypes");

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
    SnpGenotypes[] rows;
    // immutable n_ind = k_result.size1;
    // ---- Fetch annotation
    SNP[] snp_annotations;

    try {

    if (is_loco)
      snp_annotations = fetch_snp_annotations(to!string(fromStringz(cast(char *)anno_fn)));

    }
    catch (Exception e) {
      writeln(e.msg);
    }

    // ---- Parse the geno file
    info("GZipbyLine");
    info("MAF",maf_level);
    auto fn = to!string(fromStringz(cast(char *)geno_fn));
    auto use_snp_num = 0;
    bool[double] genotypes_used;
    double max_value = 0.0;
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
        if (k > max_value) max_value = k;
      }

      if (maf_level != -1) {
        // Following GEMMA logic (incorporates heterozygous information)
        // note the normalization only works when all a max value has been reached
        auto maf = elements.reduce!"a+b"/elements.length;
        if (max_value != 0.0) maf /= max_value; // normalize frequencies to between 0-1
        enforce(maf>=0.0 && maf<=1.0);
        if (maf < maf_level || maf > 1.0-maf_level) {
          trace("#snp ",line," maf:",maf," dropping ",freqs," max:",max_value);
          continue;
        }
      }
      // filter on poly
      // FIXME? The code in gemma just checks for multiple types

      // FIXME: CalcHWE

      // FIXME: r2_level

      // Scale genotypes
      auto geno_mean = reduce!"a + b"(0.0, elements)/elements.length;
      // writeln(geno_mean);
      auto elements2 = elements.map!(a => a - geno_mean);
      auto entry = (use_snp_num >= snp_annotations.length ? SnpGenotypes(NullSNP, array(elements2)) :
                    SnpGenotypes(snp_annotations[use_snp_num], array(elements2)));
      rows ~= entry;
      use_snp_num++;
    }
    info("Genotypes used: ",genotypes_used.keys.sort);
    enforce(is_centered); // FIXME
    enforce(use_snp_num > token_num,"Not enough snps to compute K "~to!string(use_snp_num));

    info("flmmd parsed ",fn," ",use_snp_num," genotypes");
    info("flmmd computes K on ",rows[0].length," individuals");

    if (is_loco) {
      auto chromosomes  = array(snp_annotations.map!(snp => snp.chr)).sort.uniq;
      info("flmmd LOCO on ",chromosomes);
      foreach(chr ; chromosomes) {
        auto outfn = to!string(fromStringz(cast(char *)target));
        outfn ~= "." ~ chr ~ (is_centered ? ".cXX.txt" : ".sXX.txt");
        auto task = task!compute_kinship(outfn,rows,chr,is_centered);
        taskPool.put(task);
      }
    }
    else {
      // ---- Compute K
      DMatrix G = new DMatrix(array(rows.map!(r => r.genotypes)));
      auto K = kinship_full(G);

      // ---- Scale K
      foreach (ref e ; K.elements) {
        e *= (1.0/G.rows);
      }
      // ---- Write K
      auto outfn = to!string(fromStringz(cast(char *)target));
      outfn ~= (is_centered ? ".cXX.txt" : ".sXX.txt");
      write_to_file(outfn,cast(immutable DMatrix)K);
    }
  }

} // C++

void compute_kinship(string outfn, SnpGenotypes[] rows, string chr, bool is_centered) {
  DMatrix G = new DMatrix(array(rows.filter!(r => r.snp.chr != chr ).map!(r => r.genotypes)));
  info("---> G",chr,":",G.rows,"x", G.cols);
  auto K = kinship_full(G);

  // ---- Scale K
  foreach (ref e ; K.elements) {
    e *= (1.0/G.rows);
  }
  // ---- Write K
  write_to_file(outfn,cast(immutable DMatrix)K);
}
