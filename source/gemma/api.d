// GEMMA C API

module gemma.api;

extern (C) {


import std.string : toStringz;


immutable(char)* faster_lmm_d_api_version() {
  return toStringz("0.1");
}

} // extern C
