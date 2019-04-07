#ifndef ZHWK_MTRACE_H
#define ZHWK_MTRACE_H

#include "flog.h"

extern FLogFile mtrace_file;

void init_mtrace(const char* mfile_path);
void term_mtrace();


#endif