#include "pthread_wrapper.h"
#include <cstdlib>
#include <errno.h>

int mu_ctor(pthread_mutex_t** mu){
    if(mu == NULL) return -EINVAL;
    (*mu) = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));
    if((*mu) == NULL) return -ENOMEM;
    return pthread_mutex_init((*mu), NULL);
}

int mu_lock(pthread_mutex_t* mu){
    if(mu == NULL) return -EINVAL;
    return pthread_mutex_lock(mu);
}

int mu_unlock(pthread_mutex_t* mu){
    if(mu == NULL) return -EINVAL;
    return pthread_mutex_unlock(mu);
}

int mu_dtor(pthread_mutex_t** mu){
    if(mu == NULL) return -EINVAL;
    free(*mu);
    (*mu) = NULL;
    return 0;
}