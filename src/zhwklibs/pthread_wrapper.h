#ifndef ZHWK_PW_H
#define ZHWK_PW_H

#include <pthread.h>


int mu_ctor(pthread_mutex_t** mu);

int mu_lock(pthread_mutex_t* mu);
int mu_unlock(pthread_mutex_t* mu);

int mu_dtor(pthread_mutex_t** mu);

#endif