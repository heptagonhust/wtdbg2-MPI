#ifndef ZHWK_FLOG_H
#define ZHWK_FLOG_H

#include <cstdio>
#include <sys/stat.h>
#include <fcntl.h>
#include <pthread.h>
#include "pthread_wrapper.h"
#include <unistd.h>

struct _flog_outfile {
    pthread_mutex_t* mu;
    FILE* fp;
};

typedef struct _flog_outfile FLogFile;

#define fLog(ffile, fmt, ...) do{ \
    mu_lock(ffile.mu);\
    fprintf(ffile.fp, fmt, __VA_ARGS__);\
    fflush(ffile.fp);
    sync();
    mu_unlock(ffile.mu);\
}while(0);

FLogFile ffile_ctor(const char* fpath);
void ffile_dtor(FLogFile* ffile);


#endif