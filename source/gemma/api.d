// GEMMA C API

module gemma.api;

import bio.std.decompress;
import std.conv;
import std.exception;
import std.stdio;

extern (C) {

  import std.string : toStringz;
  import core.stdc.string : strncpy;

  const VERSION = "0.10"c;

  // Buf should be preallocated and have at least a size of 16
  char *flmmd_api_version(char *buf) @nogc {
    strncpy(buf,VERSION.ptr,VERSION.length);
    buf[VERSION.length] = '\0';
    return buf;
  }

  void flmmd_compute_bimbam_K() {
    uint lines = 0;
    uint chars = 0;
    foreach(ubyte[] s; GzipbyLine!(ubyte[])("example/BXD_geno.txt.gz")) {
      // test file contains 7320 lines 4707218 characters
      write(cast(string)s);
      chars += s.length;
      lines += 1;
    }
    assert(chars == 4707218,"chars " ~ to!string(chars));
    assert(lines == 7321,"lines " ~ to!string(lines));
  }


} // extern C
