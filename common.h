#pragma once
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdlib.h>
char *canonicalize_file_name(const char *path);
