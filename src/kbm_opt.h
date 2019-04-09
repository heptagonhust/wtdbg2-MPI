#pragma once
#include "common.h"
#include "kbm.h"
EXTERN_C void deal_with_aux_kbm(KBMAux *aux, char *qtag, u4i qidx, BaseBank *rdseqs,
                     u8i seqoff, u4i seqlen);
EXTERN_C void split_FIXP_kmers_kbm(BaseBank *rdseqs, u8i offset, u4i length, u1i ksize,
                          u1i psize, u4i kmod, kmeroffv *rs[2]);
