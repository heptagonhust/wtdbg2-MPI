#include "flog.h"

FLogFile ffile_ctor(const char* fpath){
    FLogFile ff;
    mu_ctor(&(ff.mu));
    ff.fp = fopen(fpath, "w");
    return ff;
}

void ffile_dtor(FLogFile* ffile) {
    mu_dtor(&(ffile->mu));
    fclose(ffile->fp);
}