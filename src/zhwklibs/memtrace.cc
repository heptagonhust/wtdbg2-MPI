#include "memtrace.h"

FLogFile mtrace_file;

void init_mtrace(const char* mfile_path){
    mtrace_file = ffile_ctor(mfile_path);
}

void term_mtrace(){
    ffile_dtor(&mtrace_file);
}