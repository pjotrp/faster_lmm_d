// GEMMA C API

module gemma.api;

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

} // extern C
