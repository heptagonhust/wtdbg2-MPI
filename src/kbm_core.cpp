#include "kbm_defines.h"
#include "wtdbg.h"
#include <mpi.h>
KBMPar *init_kbmpar() {
    KBMPar *par;
    par = (KBMPar *)malloc(sizeof(KBMPar));
    par->rd_len_order = 0;
    par->use_kf = 0;
    par->min_bin_degree = 2;
    par->ksize = 0;
    par->psize = 21;
    par->kmax = 1000;
    par->kmin = 1;
    par->kmer_mod = 4 * KBM_N_HASH;
    par->ksampling = KBM_BSIZE;
    par->ktop = 0.05;
    par->strand_mask = 3;
    par->self_aln = 0;
    par->skip_contained = 1;
    par->max_bgap = 4;    // 4 * KBM_BIN_SIZE
    par->max_bvar = 4;
    par->max_bcnt = 0xFFFF;
    par->max_gap = 0.6;
    par->pgap = -7;
    par->pvar = -21;
    par->max_hit = 1000;
    par->min_aln = 2048 / KBM_BIN_SIZE;
    par->min_mat = 200;
    par->min_sim = 0.05;
    par->aln_var = 0.25;
    par->test_mode = 0;
    return par;
}

void free_kbmpar(KBMPar *par) {
    free(par);
}

KBM *init_kbm(KBMPar *par) {
    KBM *kbm;
    u4i i;
    kbm = new KBM;
    kbm->flags = 0;
    kbm->par = par;
    kbm->rdseqs = init_basebank();
    kbm->reads = init_kbmreadv(64);
    kbm->tag2idx = init_cuhash(1023);
    kbm->bins = init_kbmbinv(64);
    kbm->binmarks = init_bitvec(1024);
    //kbm->kfs = NULL;
    kbm->vec_bidxaux.reserve(64);
    for(i = 0; i < KBM_N_HASH; i++) kbm->hashs[i] = init_kbmhash(1023);
    for(i = 0; i < KBM_N_HASH; i++) kbm->kauxs[i] = init_kbmkauxv(64);
    return kbm;
}

void transfer_kbm(KBM *kbm, KBMPar *par,KBMPar *rpar,readv *reads,u4i *corr_mode,int world_rank){
    size_t size;
    u4i corrmode;
    MPI_Bcast(&kbm->flags, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (world_rank != 0) {
        kbm->par = init_kbmpar();
    }
    MPI_Bcast(kbm->par, sizeof(kbm->par), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (world_rank != 0) {
        par = init_kbmpar();
    }
    MPI_Bcast(par, sizeof(par), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (world_rank != 0) {
        rpar = init_kbmpar();
    }
    MPI_Bcast(rpar, sizeof(rpar), MPI_BYTE, 0, MPI_COMM_WORLD);


    if (world_rank == 0) {
        size = kbm->reads->size;
    }
    MPI_Bcast(&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    //brast corrmode
    // if(world_rank == 0){
    //     corrmode = *corr_mode;
    // }
    MPI_Bcast(corr_mode, 1 ,MPI_UINT32_T , 0, MPI_COMM_WORLD);

    if (world_rank != 0) {
        kbm->reads = init_kbmreadv(size);
    }
    MPI_Bcast(kbm->reads->buffer, size * sizeof(kbm_read_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    kbm->reads->size = size;
    //broadcast kbm reads

    if (world_rank == 0) {
        size = kbm->bins->size;
    }
    MPI_Bcast(&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (world_rank != 0) {
        kbm->bins = init_kbmbinv(size);
    }
    MPI_Bcast(kbm->bins->buffer, size * sizeof(kbm_bin_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    kbm->bins->size = size;
    //broadcast kbm bins

    if (world_rank == 0) {
        size = kbm->vec_bidxaux.size();
//        fprintf(stderr, "vector size: %3d: %d\n", world_rank, size);
//        fflush(stderr);
    }
    MPI_Bcast(&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    fprintf(stderr, "vector size: %3d: %d\n", world_rank, size);
    if (world_rank != 0) {
        kbm->vec_bidxaux.resize(size);
    }
    MPI_Bcast(kbm->vec_bidxaux.data(), size * sizeof(kbm_bidxaux_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    //broadcast kbm vector

    if (world_rank != 0) {
        kbm->rdseqs = (BaseBank *) malloc(sizeof(BaseBank));
    }
    MPI_Bcast(kbm->rdseqs, sizeof(BaseBank), MPI_BYTE, 0, MPI_COMM_WORLD);
    if (world_rank != 0) {
        kbm->rdseqs->bits = (u8i *) malloc(((kbm->rdseqs->size + 31) / 32 + 1) * 8);
    }
    MPI_Bcast(kbm->rdseqs->bits, (kbm->rdseqs->size + 31) / 32, MPI_UINT64_T, 0, MPI_COMM_WORLD);

    //broadcast kbm rdseqs

    for (int i = 0; i < KBM_N_HASH; i++) {
        if (world_rank != 0) {
            kbm->hashs[i] = (kbmhash *)calloc(1, sizeof(kbmhash));
        }
        MPI_Bcast(kbm->hashs[i], sizeof(kbmhash), MPI_BYTE, 0, MPI_COMM_WORLD);
        if(world_rank != 0){
            kbm->hashs[i]->array = (kbm_kmer_t *) malloc((kbm->hashs[i]->count +1) * sizeof(kbm_kmer_t));
        }
        MPI_Bcast(kbm->hashs[i]->array, kbm->hashs[i]->count + 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        //broadcast kbm hash

        if (world_rank != 0) {
            kbm->hashs[i]->ones = (BitVec *)calloc(1, sizeof(BitVec));
        }
        MPI_Bcast(kbm->hashs[i]->ones, sizeof(BitVec), MPI_BYTE, 0, MPI_COMM_WORLD);
        if (world_rank != 0){
            kbm->hashs[i]->ones->bits = (u8i*)malloc((kbm->hashs[i]->ones->n_cap / 64 + 1 )* 8);
        }
        MPI_Bcast(kbm->hashs[i]->ones->bits, kbm->hashs[i]->ones->n_cap / 64 + 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        if (world_rank != 0){
            if(kbm->hashs[i]->ones->sums)
                kbm->hashs[i]->ones->sums = (u8i*)malloc((kbm->hashs[i]->ones->sum_size * 2 + 1) * 8);
        }
        if(kbm->hashs[i]->ones->sums)
            MPI_Bcast(kbm->hashs[i]->ones->sums, kbm->hashs[i]->ones->sum_size * 2 + 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            if(kbm->hashs[i]->ones->hash)
                kbm->hashs[i]->ones->hash = (u8i*)malloc(kbm->hashs[i]->ones->hash_size * 8);
        }
        if(kbm->hashs[i]->ones->hash)
            MPI_Bcast(kbm->hashs[i]->ones->hash, kbm->hashs[i]->ones->hash_size, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        //broadcast kbm hash one

        if(kbm->hashs[i]->dels) {
            if (world_rank != 0) {
                kbm->hashs[i]->dels = (BitVec *)calloc(1, sizeof(BitVec));
            }
            MPI_Bcast(kbm->hashs[i]->dels, sizeof(BitVec), MPI_BYTE, 0, MPI_COMM_WORLD);
            if (world_rank != 0){
                kbm->hashs[i]->dels->bits = (u8i*)malloc((kbm->hashs[i]->dels->n_cap / 64 + 1 )* 8);
            }
            MPI_Bcast(kbm->hashs[i]->dels->bits, kbm->hashs[i]->dels->n_cap / 64 + 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
            if (world_rank != 0){
                if(kbm->hashs[i]->dels->sums)
                    kbm->hashs[i]->dels->sums = (u8i*)malloc((kbm->hashs[i]->dels->sum_size * 2 + 1) * 8);
            }
            if(kbm->hashs[i]->dels->sums)
                MPI_Bcast(kbm->hashs[i]->dels->sums, kbm->hashs[i]->dels->sum_size * 2 + 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
            if (world_rank != 0) {
                if(kbm->hashs[i]->dels->hash)
                    kbm->hashs[i]->dels->hash = (u8i*)malloc(kbm->hashs[i]->dels->hash_size * 8);
            }
            if(kbm->hashs[i]->dels->hash)
                MPI_Bcast(kbm->hashs[i]->dels->hash, kbm->hashs[i]->dels->hash_size, MPI_UINT64_T, 0, MPI_COMM_WORLD);
            //broadcast kbm hash del
        }

        if (world_rank == 0) {
            size = kbm->kauxs[i]->size;
        }
        MPI_Bcast(&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            kbm->kauxs[i] = init_kbmkauxv(size);
        }
        MPI_Bcast(kbm->kauxs[i]->buffer, size, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        kbm->kauxs[i]->size = size;
        //broadcast kbm kauxs
    }
    //if (corr_mode) trans g->reads
    if(*corr_mode){
         if(world_rank != 0){
             reads = init_readv(kbm->reads->size);
         }
         MPI_Bcast(reads , kbm->reads->size*sizeof(read_t), MPI_BYTE, 0 ,MPI_COMM_WORLD);
    }

}

void free_kbm(KBM *kbm) {
    u4i i;
    if(kbm->flags & 0x1) return;
    if(kbm->flags & 0x2) {
    } else {
        free_basebank(kbm->rdseqs);
        for(i = 0; i < kbm->reads->size; i++) {
            if(kbm->reads->buffer[i].tag) free(kbm->reads->buffer[i].tag);
        }
        free_kbmreadv(kbm->reads);
    }
    free_cuhash(kbm->tag2idx);
    free_kbmbinv(kbm->bins);
    free_bitvec(kbm->binmarks);
    //if(kbm->kfs) free(kbm->kfs);
    for(i = 0; i < KBM_N_HASH; i++) free_kbmhash(kbm->hashs[i]);
    for(i = 0; i < KBM_N_HASH; i++) free_kbmkauxv(kbm->kauxs[i]);
    delete kbm;
}

void reset_index_kbm(KBM *kbm) {
    u4i i;
    if(kbm->flags & 0x04) {
        free_kbmhash(kbm->hashs[0]);
        kbm->hashs[0] = init_kbmhash(1023);
        free_kbmkauxv(kbm->kauxs[0]);
        kbm->kauxs[0] = init_kbmkauxv(64);
    } else {
        for(i = 0; i < KBM_N_HASH; i++) {
            free_kbmhash(kbm->hashs[i]);
            kbm->hashs[i] = init_kbmhash(1023);
            free_kbmkauxv(kbm->kauxs[i]);
            kbm->kauxs[i] = init_kbmkauxv(64);
        }
    }
    // free_kbmbmerv(kbm->seeds);
    // kbm->seeds = init_kbmbmerv(64);
    // free_kbmbauxv(kbm->sauxs);
    // kbm->sauxs = init_kbmbauxv(64);
    kbm->vec_bidxaux.clear();
    //zeros_bitvec(kbm->binmarks);
}

void clear_kbm(KBM *kbm) {
    u4i i;
    if(kbm->flags & 0x1) return;
    if(kbm->flags & 0x2) {
    } else {
        clear_basebank(kbm->rdseqs);
        for(i = 0; i < kbm->reads->size; i++) {
            if(kbm->reads->buffer[i].tag) free(kbm->reads->buffer[i].tag);
        }
        clear_kbmreadv(kbm->reads);
        clear_cuhash(kbm->tag2idx);
    }
    clear_kbmbinv(kbm->bins);
    clear_bitvec(kbm->binmarks);
    for(i = 0; i < KBM_N_HASH; i++) {
        free_kbmhash(kbm->hashs[i]);
        kbm->hashs[i] = init_kbmhash(1023);
        free_kbmkauxv(kbm->kauxs[i]);
        kbm->kauxs[i] = init_kbmkauxv(64);
        if(kbm->flags & 0x04) break;
    }
    // clear_kbmbmerv(kbm->seeds);
    // clear_kbmbauxv(kbm->sauxs);
    kbm->vec_bidxaux.clear();
}

int cvt_kbm_read_length(u4i seqlen) {
    if(seqlen > KBM_MAX_RDLEN) seqlen = KBM_MAX_RDLEN;
    seqlen = (seqlen / KBM_BIN_SIZE) * KBM_BIN_SIZE;
    return seqlen;
}

void push_kbm(KBM *kbm, char *tag, int taglen, char *seq, u4i seqlen) {
    kbm_read_t *rd;
    char *ptr;
    if(taglen) {
        ptr = (char *)malloc(taglen + 1);
        memcpy(ptr, tag, taglen);
        ptr[taglen] = 0;
    } else {
        ptr = NULL;
    }
    seqlen = cvt_kbm_read_length(seqlen);
    rd = next_ref_kbmreadv(kbm->reads);
    rd->rdoff = kbm->rdseqs->size;
    rd->rdlen = seqlen;
    rd->binoff = 0;
    rd->bincnt = 0;
    rd->tag = ptr;
    seq2basebank(kbm->rdseqs, seq, seqlen);
    // make sure rdoff is even
    if(kbm->rdseqs->size & 0x1) {
        bit2basebank(kbm->rdseqs, 0);
    }
}

void bitpush_kbm(KBM *kbm, char *tag, int taglen, u8i *seqs, u8i seqoff, u4i seqlen) {
    kbm_read_t *rd;
    char *ptr;
    if(taglen) {
        ptr = (char *)malloc(taglen + 1);
        memcpy(ptr, tag, taglen);
        ptr[taglen] = 0;
    } else {
        ptr = NULL;
    }
    seqlen = cvt_kbm_read_length(seqlen);
    rd = next_ref_kbmreadv(kbm->reads);
    rd->rdoff = kbm->rdseqs->size;
    rd->rdlen = seqlen;
    rd->binoff = 0;
    rd->bincnt = 0;
    rd->tag = ptr;
    fast_fwdbits2basebank(kbm->rdseqs, seqs, seqoff, seqlen);
    // make sure rdoff is even
    if(kbm->rdseqs->size & 0x1) {
        bit2basebank(kbm->rdseqs, 0);
    }
}

// Please call no more than once
u8i filter_reads_kbm(KBM *kbm, u8i retain_size, int strategy) {
    u8i m, b, e, len;
    if(kbm->reads->size == 0) return 0;
    if(retain_size == 0 || retain_size >= kbm->rdseqs->size) return kbm->rdseqs->size;
    if((kbm->flags & 0x2) == 0) {
        if(kbm->par->rd_len_order) {
            sort_array(kbm->reads->buffer, kbm->reads->size, kbm_read_t,
                       num_cmpgt(b.rdlen, a.rdlen));
            if(strategy == 0) {    // longest
                len = 0;
                for(e = 0; e < kbm->reads->size; e++) {
                    len += kbm->reads->buffer[e].rdlen;
                    if(len >= retain_size) break;
                }
                kbm->reads->size = e;
            } else if(strategy == 1) {    // median
                m = kbm->reads->size / 2;
                len = kbm->reads->buffer[m].rdlen;
                e = m;
                for(b = 0; b <= m && len < retain_size; b++) {
                    len += kbm->reads->buffer[m - b].rdlen;
                    len += kbm->reads->buffer[m + b].rdlen;
                }
                e = b * 2;
                b = m - b;
                if(b) {
                    remove_array_kbmreadv(kbm->reads, 0, b);
                }
                kbm->reads->size = e;
            } else {
                return kbm->rdseqs->size;
            }
            return len;
        } else {
            return kbm->rdseqs->size;
        }
    } else {
        return kbm->rdseqs->size;
    }
}

void ready_kbm(KBM *kbm) {
    kbm_read_t *rd;
    u4i i, j;
    if((kbm->flags & 0x2) == 0) {
        if(kbm->par->rd_len_order) {
            sort_array(kbm->reads->buffer, kbm->reads->size, kbm_read_t,
                       num_cmpgt(b.rdlen, a.rdlen));
        }
        encap_basebank(kbm->rdseqs, KBM_BSIZE);
    }
    clear_kbmbinv(kbm->bins);
    for(i = 0; i < kbm->reads->size; i++) {
        rd = ref_kbmreadv(kbm->reads, i);
        if(rd->tag) put_cuhash(kbm->tag2idx, (cuhash_t){rd->tag, i});
        if((kbm->flags & 0x2) == 0) rd->binoff = kbm->bins->size;
        for(j = 0; j + KBM_BIN_SIZE <= rd->rdlen; j += KBM_BIN_SIZE) {
            push_kbmbinv(kbm->bins, (kbm_bin_t){i, j / KBM_BIN_SIZE, 0, 0});
        }
        if((kbm->flags & 0x2) == 0) rd->bincnt = j / KBM_BIN_SIZE;
    }
    clear_bitvec(kbm->binmarks);
    encap_bitvec(kbm->binmarks, kbm->bins->size);
    kbm->binmarks->n_bit = kbm->bins->size;
    zeros_bitvec(kbm->binmarks);
}

// Share seqs, reads
KBM *clone_seqs_kbm(KBM *src, KBMPar *par) {
    KBM *dst;
    dst = init_kbm(par);
    free_basebank(dst->rdseqs);
    free_kbmreadv(dst->reads);
    dst->rdseqs = src->rdseqs;
    dst->reads = src->reads;
    dst->flags = 1LLU << 1;
    ready_kbm(dst);    // Notice encap_basebank in ready_kbm
    return dst;
}


u8i seed2solid_idx_kbm(KBM *kbm, kbm_dpe_t *p) {
    kbm_bin_t *b;
    kbm_read_t *rd;
    u8i seqoff;
    b = kbm->bins->buffer + p->bidx;
    rd = kbm->reads->buffer + b->ridx;
    seqoff = ((rd->rdoff + b->off * KBM_BSIZE) >> 1) + p->koff;
    return seqoff;
}

u8i rdoff2solid_idx_kbm(KBM *kbm, u4i ridx, u4i roff) {
    kbm_read_t *rd;
    u8i seqoff;
    rd = kbm->reads->buffer + ridx;
    seqoff = (rd->rdoff + roff) >> 1;
    return seqoff;
}

#define binoff2solid_koff_kbm(kbm, bidx, boff) ((boff) >> 1)

typedef struct {
    u8i mer : 50, kidx : 14;
    u8i bidx;
    u4i cnt : 22, koff : 8, dir : 1, used : 1;
} kbm_midx_t;
define_list(kbmmidxv, kbm_midx_t);

// typedef struct {
//     u8i bidx;
//     kbm_baux_t aux;
// } kbm_tmp_bmer_t;
// define_list(tmpbmerv, kbm_tmp_bmer_t);

thread_beg_def(midx);
KBM *kbm;
u4i beg, end;    // (end - beg) * KBM_BSIZE MUST <= KBM_KMEROFF_MAX
u8i ktot, nrem, Nrem, none, nflt, offset;
u8i srem, Srem;
u8i *cnts, n_cnt;
int task;
int cal_degree;
pthread_mutex_t *locks;
thread_end_def(midx);

thread_beg_func(midx);
KBM *kbm;
kbm_bin_t *bin;
kbmmidxv **kidxs;
kbm_midx_t *mx;
kmer_off_t *f;
kbm_kmer_t *u;
kbm_kaux_t *x;
kmeroffv *kmers[2];
u8i off;
u4i bidx, i, j, k, len, ncpu, tidx, kidx;
int exists;
kbm = midx->kbm;
ncpu = midx->n_cpu;
tidx = midx->t_idx;
kmers[0] = adv_init_kmeroffv(64, 0, 1);
kmers[1] = adv_init_kmeroffv(64, 0, 1);
kidxs = (kbmmidxv **)malloc(KBM_N_HASH * sizeof(kbmmidxv *));
for(i = 0; i < KBM_N_HASH; i++) kidxs[i] = init_kbmmidxv(64);
thread_beg_loop(midx);
if(midx->task == 1) {
    // counting kmers
    for(i = 0; i < KBM_N_HASH; i++) clear_kbmmidxv(kidxs[i]);
    for(bidx = midx->beg + tidx; bidx < midx->end; bidx += ncpu) {
        if(KBM_LOG == 0 && tidx == 0 && ((bidx - midx->beg) % 100000) == 0) {
            fprintf(KBM_LOGF, "\r%u", bidx - midx->beg);
            fflush(KBM_LOGF);
        }
        bin = ref_kbmbinv(kbm->bins, bidx);
        if(bin->closed) continue;
        if(((u4i)bin->off + 1) * KBM_BIN_SIZE > kbm->reads->buffer[bin->ridx].rdlen)
            continue;
        off = kbm->reads->buffer[bin->ridx].rdoff + bin->off * KBM_BIN_SIZE;
        len = KBM_BIN_SIZE;
        split_FIXP_kmers_kbm(kbm->rdseqs, off, len, kbm->par->ksize, kbm->par->psize,
                             kbm->par->kmer_mod, kmers);
        for(i = 0; i < 2; i++) {
            for(j = 0; j < kmers[i]->size; j++) {
                f = ref_kmeroffv(kmers[i], j);
                if(f->closed) continue;
                mx = next_ref_kbmmidxv(kidxs[f->kidx]);
                mx->mer = f->kmer;
                mx->kidx = f->kidx;
                mx->bidx = bidx;
                mx->dir = i;
                mx->koff = f->off;
                if(kidxs[f->kidx]->size >= 64) {
                    kidx = f->kidx;
                    // lock hashs[kidx]
                    pthread_mutex_lock(midx->locks + kidx);
                    // hash adding
                    for(k = 0; k < kidxs[kidx]->size; k++) {
                        mx = ref_kbmmidxv(kidxs[kidx], k);
                        u = prepare_kbmhash(kbm->hashs[kidx], mx->mer, &exists);
                        if(exists) {
                            if(u->tot < KBM_MAX_KCNT) u->tot++;
                        } else {
                            u->mer = mx->mer;
                            u->tot = 1;
                            u->flt = 0;
                        }
                    }
                    // free hashs[f->kidx]
                    pthread_mutex_unlock(midx->locks + kidx);
                    clear_kbmmidxv(kidxs[kidx]);
                }
            }
        }
    }
    for(kidx = 0; kidx < KBM_N_HASH; kidx++) {
        if(kidxs[kidx]->size) {
            // lock hashs[kidx]
            pthread_mutex_lock(midx->locks + kidx);
            // hash adding
            for(k = 0; k < kidxs[kidx]->size; k++) {
                mx = ref_kbmmidxv(kidxs[kidx], k);
                u = prepare_kbmhash(kbm->hashs[kidx], mx->mer, &exists);
                if(exists) {
                    if(u->tot < KBM_MAX_KCNT) u->tot++;
                } else {
                    u->mer = mx->mer;
                    u->tot = 1;
                    u->flt = 0;
                }
            }
            // free hashs[f->kidx]
            pthread_mutex_unlock(midx->locks + kidx);
            clear_kbmmidxv(kidxs[kidx]);
        }
    }
} else if(midx->task == 2) {
    // delete low freq kmers
    for(i = tidx; i < KBM_N_HASH; i += ncpu) {
        reset_iter_kbmhash(kbm->hashs[i]);
        while((u = ref_iter_kbmhash(kbm->hashs[i]))) {
            if(u->tot < kbm->par->kmin) {
                delete_kbmhash(kbm->hashs[i], u);
            }
        }
    }
} else if(midx->task == 3) {
    // stat kmer counts
    memset(midx->cnts, 0, midx->n_cnt * sizeof(u8i));
    for(i = tidx; i < KBM_N_HASH; i += ncpu) {
        reset_iter_kbmhash(kbm->hashs[i]);
        while((u = ref_iter_kbmhash(kbm->hashs[i]))) {
            midx->cnts[num_min(midx->n_cnt, u->tot) - 1]++;
        }
    }
} else if(midx->task == 4) {
    // stat counts
    midx->offset = 0;
    for(i = tidx; i < KBM_N_HASH; i += ncpu) {
        reset_iter_kbmhash(kbm->hashs[i]);
        while((u = ref_iter_kbmhash(kbm->hashs[i]))) {
            x = ref_kbmkauxv(kbm->kauxs[i], offset_kbmhash(kbm->hashs[i], u));
            x->off = midx->offset;
            x->cnt = 0;
            midx->ktot += u->tot;
            if(u->tot < kbm->par->kmin) {
                u->flt = 1;
            } else if(u->tot > kbm->par->kmax) {
                u->flt = 1;
            } else {
                midx->offset += u->tot;
            }
        }
    }
} else if(midx->task == 5) {
    // revise offset
    for(i = tidx; i < KBM_N_HASH; i += ncpu) {
        for(off = 0; off < kbm->kauxs[i]->size; off++)
            kbm->kauxs[i]->buffer[off].off += midx->offset;
    }
} else if(midx->task == 6) {
    // fill seeds
    for(i = 0; i < KBM_N_HASH; i++) clear_kbmmidxv(kidxs[i]);
    u4v *chgs;
    chgs = init_u4v(KBM_BSIZE);
    for(bidx = midx->beg + tidx; bidx < midx->end; bidx += ncpu) {
        if(KBM_LOG == 0 && tidx == 0 && ((bidx - midx->beg) % 100000) == 0) {
            fprintf(KBM_LOGF, "\r%u", bidx - midx->beg);
            fflush(KBM_LOGF);
        }
        bin = ref_kbmbinv(kbm->bins, bidx);
        if(bin->closed) continue;
        bin->degree = 0;
        if(((u4i)bin->off + 1) * KBM_BIN_SIZE > kbm->reads->buffer[bin->ridx].rdlen)
            continue;
        off = kbm->reads->buffer[bin->ridx].rdoff + bin->off * KBM_BIN_SIZE;
        len = KBM_BIN_SIZE;
        split_FIXP_kmers_kbm(kbm->rdseqs, off, len, kbm->par->ksize, kbm->par->psize,
                             kbm->par->kmer_mod, kmers);
        clear_u4v(chgs);
        for(i = 0; i < 2; i++) {
            for(j = 0; j < kmers[i]->size; j++) {
                f = ref_kmeroffv(kmers[i], j);
                if(f->closed) continue;
                mx = next_ref_kbmmidxv(kidxs[f->kidx]);
                mx->mer = f->kmer;
                mx->kidx = f->kidx;
                mx->bidx = bidx;
                mx->dir = i;
                mx->koff = f->off;
                if(kidxs[f->kidx]->size == 64) {
                    push_u4v(chgs, f->kidx);
                }
            }
        }
        for(i = 0; i < chgs->size; i++) {
            kidx = chgs->buffer[i];
            if(kidxs[kidx]->size == 0) continue;
            pthread_mutex_lock(midx->locks + kidx);
            for(k = 0; k < kidxs[kidx]->size; k++) {
                mx = ref_kbmmidxv(kidxs[kidx], k);
                u = get_kbmhash(kbm->hashs[kidx], mx->mer);
                if(u && u->flt == 0) {
                    x = ref_kbmkauxv(kbm->kauxs[kidx],
                                     offset_kbmhash(kbm->hashs[kidx], u));
                    kbm->bins->buffer[mx->bidx].degree++;
                    if(x->cnt < u->tot) {
                        if(x->cnt && getval_bidx(kbm, x->off + x->cnt - 1) == mx->bidx &&
                           kbm->vec_bidxaux[x->off + x->cnt - 1].dir == mx->dir) {
                            // repeated kmer within one bin
                        } else {
                            kbm->vec_bidxaux[x->off + x->cnt].bidx = mx->bidx;
                            kbm->vec_bidxaux[x->off + x->cnt].dir = mx->dir;
                            kbm->vec_bidxaux[x->off + x->cnt].koff = mx->koff >> 1;
                            x->cnt++;
                        }
                    }
                }
            }
            pthread_mutex_unlock(midx->locks + kidx);
            clear_kbmmidxv(kidxs[kidx]);
        }
    }
    for(kidx = 0; kidx < KBM_N_HASH; kidx++) {
        if(kidxs[kidx]->size) {
            // lock hashs[kidx]
            pthread_mutex_lock(midx->locks + kidx);
            // hash adding
            for(k = 0; k < kidxs[kidx]->size; k++) {
                mx = ref_kbmmidxv(kidxs[kidx], k);
                u = get_kbmhash(kbm->hashs[kidx], mx->mer);
                if(u && u->flt == 0) {
                    x = ref_kbmkauxv(kbm->kauxs[kidx],
                                     offset_kbmhash(kbm->hashs[kidx], u));
                    kbm->bins->buffer[mx->bidx].degree++;
                    if(x->cnt < u->tot) {
                        if(x->cnt && getval_bidx(kbm, x->off + x->cnt - 1) == mx->bidx &&
                           kbm->vec_bidxaux[x->off + x->cnt - 1].dir == mx->dir) {
                            // repeated kmer within one bin
                        } else {
                            kbm->vec_bidxaux[x->off + x->cnt].bidx = mx->bidx;
                            kbm->vec_bidxaux[x->off + x->cnt].dir = mx->dir;
                            kbm->vec_bidxaux[x->off + x->cnt].koff = mx->koff >> 1;
                            x->cnt++;
                        }
                    }
                }
            }
            // free hashs[f->kidx]
            pthread_mutex_unlock(midx->locks + kidx);
            clear_kbmmidxv(kidxs[kidx]);
        }
    }
    free_u4v(chgs);
} else if(midx->task == 7) {
    // count added kmers
    midx->srem = midx->Srem = 0;
    for(i = tidx; i < KBM_N_HASH; i += ncpu) {
        reset_iter_kbmhash(kbm->hashs[i]);
        while((u = ref_iter_kbmhash(kbm->hashs[i]))) {
            x = ref_kbmkauxv(kbm->kauxs[i], offset_kbmhash(kbm->hashs[i], u));
            if(x->cnt) {
                midx->srem++;
                midx->Srem += x->cnt;
            }
        }
    }
} else if(midx->task == 8) {
    // sort seeds within a kmer
    for(i = tidx; i < KBM_N_HASH; i += ncpu) {
        reset_iter_kbmhash(kbm->hashs[i]);
        while((u = ref_iter_kbmhash(kbm->hashs[i]))) {
            x = ref_kbmkauxv(kbm->kauxs[i], offset_kbmhash(kbm->hashs[i], u));
            if(x->cnt < 2) continue;
            auto beg_iter = kbm->vec_bidxaux.begin() + x->off;
            auto end_iter = beg_iter + x->cnt;
            std::sort(beg_iter, end_iter, [](auto& a, auto& b){return a.bidx < b.bidx;});
            // for(j = 0; j < x->cnt; j++) {
            //     push_tmpbmerv(bms, (kbm_tmp_bmer_t){getval_bidx(kbm, x->off + j),
            //                                         kbm->sauxs->buffer[x->off + j]});
            // }
            // sort_array(bms->buffer, bms->size, kbm_tmp_bmer_t, num_cmpgt(a.bidx, b.bidx));
            // kbm->seeds->buffer[x->off + 0].bidx = bms->buffer[0].bidx & MAX_U4;
            // kbm->sauxs->buffer[x->off + 0].bidx = bms->buffer[0].bidx >> 32;
            // kbm->sauxs->buffer[x->off + 0] = bms->buffer[0].aux;
            // len = 1;
            // for(j = 1; j < x->cnt; j++) {
            //     // if(bms->buffer[j].bidx < bms->buffer[j - 1].bidx) {
            //     //     continue;
            //     // }
            //     kbm->seeds->buffer[x->off + len].bidx = bms->buffer[j].bidx & MAX_U4;
            //     kbm->sauxs->buffer[x->off + len].bidx = bms->buffer[0].bidx >> 32;
            //     kbm->sauxs->buffer[x->off + len] = bms->buffer[j].aux;
            //     len++;
            // }
            // x->cnt = len;
        }
    }
}
thread_end_loop(midx);
free_kmeroffv(kmers[0]);
free_kmeroffv(kmers[1]);
for(i = 0; i < KBM_N_HASH; i++) free_kbmmidxv(kidxs[i]);
free(kidxs);
thread_end_func(midx);

void index_kbm(KBM *kbm, u8i beg, u8i end, u4i ncpu, FILE *kmstat) {
    u8i ktyp, nflt, nrem, Nrem, none, ktot, srem, Srem, off, cnt, *kcnts, MAX;
    u8i i, b, e, n;
    u4i j, kavg, batch_size;
    pthread_mutex_t *hash_locks;
    thread_preprocess(midx);
    batch_size = 10000;
    // DANGEROUS
    // clear_kbmbmerv(kbm->seeds);
    kbm->vec_bidxaux.clear();
    for(i = 0; i < KBM_N_HASH; i++) clear_kbmhash(kbm->hashs[i]);
    kcnts = NULL;
    MAX = KBM_MAX_KCNT;
    //if(kbm->kfs){
    //free(kbm->kfs);
    //kbm->kfs = NULL;
    //}
    //if(kbm->par->kmin <= 1) kbm->par->use_kf = 0;
    //kbm->kfs = kbm->par->use_kf? calloc(KBM_KF_SIZE / 4 / 8, 8) : NULL;
    hash_locks = (pthread_mutex_t *)calloc(KBM_N_HASH, sizeof(pthread_mutex_t));
    thread_beg_init(midx, ncpu);
    midx->kbm = kbm;
    midx->beg = beg;
    midx->end = end;
    midx->cnts = NULL;
    midx->n_cnt = MAX;
    midx->task = 0;
    midx->cal_degree = 0;
    midx->locks = hash_locks;
    thread_end_init(midx);
    fprintf(KBM_LOGF, "[%s] - scanning kmers (K%dP%dS%0.2f) from %llu bins\n", date(),
            kbm->par->ksize, kbm->par->psize, 1.0 * kbm->par->kmer_mod / KBM_N_HASH,
            end - beg);
    b = e = beg;
    thread_apply_all(midx, midx->task = 1);
    if(KBM_LOG == 0) {
        fprintf(KBM_LOGF, "\r%llu bins\n", end - beg);
        fflush(KBM_LOGF);
    }
    kcnts = (u8i *)calloc(MAX, sizeof(u8i));
    thread_beg_iter(midx);
    midx->cnts = (u8i *)calloc(MAX, sizeof(u8i));
    midx->task = 3;    // counting raw kmers
    thread_wake(midx);
    thread_end_iter(midx);
    thread_beg_iter(midx);
    thread_wait(midx);
    for(i = 0; i < MAX; i++) kcnts[i] += midx->cnts[i];
    thread_end_iter(midx);
    if(kmstat) {
        fprintf(kmstat, "#Reads: %llu\n", (u8i)kbm->reads->size);
        fprintf(kmstat, "#Bases: %llu bp\n", (u8i)kbm->rdseqs->size);
        fprintf(kmstat, "#K%dP%dS%0.2f\n", kbm->par->ksize, kbm->par->psize,
                1.0 * kbm->par->kmer_mod / KBM_N_HASH);
        for(i = 0; i + 1 < MAX; i++) {
            fprintf(kmstat, "%llu\t%llu\t%llu\n", i + 1, kcnts[i], (i + 1) * kcnts[i]);
        }
        fflush(kmstat);
    }
    if(kmstat) {
        u8i *_kcnts;
        _kcnts = (u8i *)malloc(200 * sizeof(u8i));
        for(i = 0; i < 200; i++) {
            _kcnts[i] = (i + 1) * kcnts[i];
        }
        char *txt = barplot_txt_u8_simple(100, 20, _kcnts, 200, 0);
        fprintf(KBM_LOGF,
                "********************** Kmer Frequency **********************\n");
        fputs(txt, KBM_LOGF);
        fprintf(KBM_LOGF,
                "**********************     1 - 201    **********************\n");
        ktyp = 0;
        for(i = 0; i < MAX; i++) {
            ktyp += (i + 1) * kcnts[i];
        }
        _kcnts[0] = 0.10 * ktyp;
        _kcnts[1] = 0.20 * ktyp;
        _kcnts[2] = 0.30 * ktyp;
        _kcnts[3] = 0.40 * ktyp;
        _kcnts[4] = 0.50 * ktyp;
        _kcnts[5] = 0.60 * ktyp;
        _kcnts[6] = 0.70 * ktyp;
        _kcnts[7] = 0.80 * ktyp;
        _kcnts[8] = 0.90 * ktyp;
        _kcnts[9] = 0.95 * ktyp;
        fprintf(KBM_LOGF, "Quatiles:\n");
        fprintf(
            KBM_LOGF,
            "   10%%   20%%   30%%   40%%   50%%   60%%   70%%   80%%   90%%   95%%\n");
        off = 0;
        for(i = j = 0; j < 10; j++) {
            while(off < _kcnts[j] && i < MAX) {
                off += kcnts[i] * (i + 1);
                i++;
            }
            fprintf(KBM_LOGF, "%6llu", i);
        }
        fprintf(KBM_LOGF, "\n");
        fprintf(KBM_LOGF,
                "# If the kmer distribution is not good, please kill me and adjust -k, "
                "-p, and -K\n"
                "# Cannot get a good distribution anyway, should adjust -S -s, also -A "
                "-e in assembly\n");
        free(_kcnts);
        free(txt);
    }
    // delete low freq kmer from hash
    thread_apply_all(midx, midx->task = 2);
    // freeze hash to save memory and speed up the query
    for(i = 0; i < KBM_N_HASH; i++) {
        if(0) {
            fprintf(KBM_LOGF, "%12llu ", (u8i)kbm->hashs[i]->count);
            fflush(KBM_LOGF);
            if((i % 8) == 7) {
                fprintf(KBM_LOGF, "\n");
                fflush(KBM_LOGF);
            }
        }
        freeze_kbmhash(kbm->hashs[i], 1.0 / 16);
        free_kbmkauxv(kbm->kauxs[i]);
        kbm->kauxs[i] = init_kbmkauxv(kbm->hashs[i]->count);
        kbm->kauxs[i]->size = kbm->hashs[i]->count;
    }
    print_proc_stat_info(0);
    ktot = nrem = Nrem = none = nflt = ktyp = 0;
    for(i = 0; i + 1 < MAX; i++) {
        ktot += kcnts[i];
        if(i + 1 < kbm->par->kmin) {
            none += kcnts[i];
        } else {
            ktyp += kcnts[i] * (i + 1);
        }
    }
    if(kbm->par->ktop) {
        off = 0;
        for(i = MAX; i > kbm->par->kmin; i--) {
            off += kcnts[i - 1] * i;
            if(off >= (ktyp * kbm->par->ktop)) break;
        }
        if(i > kbm->par->kmax) {
            kbm->par->kmax = i;
        }
        fprintf(KBM_LOGF, "[%s] - high frequency kmer depth is set to %d\n", date(),
                kbm->par->kmax);
    }
    for(i = kbm->par->kmin ? kbm->par->kmin - 1 : 0; i < kbm->par->kmax; i++) {
        nrem += kcnts[i];
        Nrem += kcnts[i] * (i + 1);
    }
    for(i = kbm->par->kmax; i + 1 < MAX; i++) {
        nflt += kcnts[i];
    }
    off = 0;
    thread_apply_all(midx, midx->task = 4);
    thread_beg_iter(midx);
    cnt = midx->offset;
    midx->offset = off;
    off += cnt;
    midx->task = 5;
    thread_wake(midx);
    thread_end_iter(midx);
    thread_wait_all(midx);
    // clear_and_encap_kbmbmerv(kbm->seeds, off + 1);
    // kbm->seeds->size = off;
    // free_kbmbauxv(kbm->sauxs);
    // kbm->sauxs = init_kbmbauxv(off + 1);
    // kbm->sauxs->size = off;
    kbm->vec_bidxaux.resize(off);
    kavg = Nrem / (nrem + 1);
    fprintf(KBM_LOGF, "[%s] - Total kmers = %llu\n", date(), ktot);
    fprintf(KBM_LOGF, "[%s] - average kmer depth = %d\n", date(), kavg);
    fprintf(KBM_LOGF, "[%s] - %llu low frequency kmers (<%d)\n", date(), none,
            kbm->par->kmin);
    fprintf(KBM_LOGF, "[%s] - %llu high frequency kmers (>%d)\n", date(), nflt,
            kbm->par->kmax);
    fprintf(KBM_LOGF, "[%s] - indexing %llu kmers, %llu instances (at most)\n", date(),
            nrem, Nrem);
    thread_apply_all(midx, midx->task = 6);
    if(KBM_LOG == 0) {
        fprintf(KBM_LOGF, "\r%llu bins\n", end - beg);
        fflush(KBM_LOGF);
    }
    thread_apply_all(midx, midx->task = 7);
    srem = Srem = 0;
    thread_beg_iter(midx);
    srem += midx->srem;
    Srem += midx->Srem;
    thread_end_iter(midx);
    fprintf(KBM_LOGF, "[%s] - indexed  %llu kmers, %llu instances\n", date(), srem, Srem);
    {
        n = 0;
        for(i = beg; i < end; i++) {
            if(kbm->bins->buffer[i].degree < kbm->par->min_bin_degree) {
                kbm->bins->buffer[i].closed = 1;
                one_bitvec(kbm->binmarks, i);
                n++;
            }
        }
        index_bitvec(kbm->binmarks);
        fprintf(KBM_LOGF, "[%s] - masked %llu bins as closed\n", date(), n);
    }
    fprintf(KBM_LOGF, "[%s] - sorting\n", date());
    thread_apply_all(midx, midx->task = 8);
    if(kcnts) free(kcnts);
    thread_beg_close(midx);
    if(midx->cnts) free(midx->cnts);
    thread_end_close(midx);
    if(1) {
        free(hash_locks);
    }
    print_proc_stat_info(0);
}

void simple_index_kbm(KBM *kbm, u8i beg, u8i end) {
    kbm_bin_t *bin;
    kbm_kmer_t *u;
    kbm_kaux_t *x;
    kmeroffv *kmers[2];
    // tmpbmerv *bms;
    kmer_off_t *f;
    u8i bidx, off;
    u4i i, j, len;
    int exists;
    kbm->flags |= 1LLU << 2;
    kmers[0] = adv_init_kmeroffv(64, 0, 1);
    kmers[1] = adv_init_kmeroffv(64, 0, 1);
    // bms = init_tmpbmerv(KBM_MAX_KCNT);
    // count kmers
    {
        for(bidx = beg; bidx < end; bidx++) {
            bin = ref_kbmbinv(kbm->bins, bidx);
            if(bin->closed) continue;
            if(((u4i)bin->off + 1) * KBM_BIN_SIZE > kbm->reads->buffer[bin->ridx].rdlen)
                continue;
            off = kbm->reads->buffer[bin->ridx].rdoff + bin->off * KBM_BIN_SIZE;
            split_FIXP_kmers_kbm(kbm->rdseqs, off, KBM_BIN_SIZE, kbm->par->ksize,
                                 kbm->par->psize, kbm->par->kmer_mod, kmers);
            for(i = 0; i < 2; i++) {
                for(j = 0; j < kmers[i]->size; j++) {
                    f = ref_kmeroffv(kmers[i], j);
                    if(f->closed) continue;
                    u = prepare_kbmhash(kbm->hashs[0], f->kmer, &exists);
                    if(exists) {
                        if(u->tot < KBM_MAX_KCNT) u->tot++;
                    } else {
                        u->mer = f->kmer;
                        u->tot = 1;
                        u->flt = 0;
                    }
                }
            }
        }
    }
    // delete low freq kmers
    if(kbm->par->kmin) {
        reset_iter_kbmhash(kbm->hashs[0]);
        while((u = ref_iter_kbmhash(kbm->hashs[0]))) {
            if(u->tot < kbm->par->kmin) {
                delete_kbmhash(kbm->hashs[0], u);
            }
        }
    }
    // freeze hash
    {
        freeze_kbmhash(kbm->hashs[0], 1.0 / 16);
        free_kbmkauxv(kbm->kauxs[0]);
        kbm->kauxs[0] = init_kbmkauxv(kbm->hashs[0]->count);
        kbm->kauxs[0]->size = kbm->hashs[0]->count;
    }
    // encap seeds
    {
        off = 0;
        reset_iter_kbmhash(kbm->hashs[0]);
        while((u = ref_iter_kbmhash(kbm->hashs[0]))) {
            x = ref_kbmkauxv(kbm->kauxs[0], offset_kbmhash(kbm->hashs[0], u));
            x->off = off;
            x->cnt = 0;
            if(u->tot < kbm->par->kmin) {    // Never happens
                u->flt = 1;
            } else if(kbm->par->kmax && u->tot > kbm->par->kmax) {
                u->flt = 1;
            } else {
                off += u->tot;
            }
        }
        // clear_and_encap_kbmbmerv(kbm->seeds, off + 1);
        // kbm->seeds->size = off;
        // clear_and_encap_kbmbauxv(kbm->sauxs, off + 1);
        // kbm->sauxs->size = off;
        kbm->vec_bidxaux.clear();
        kbm->vec_bidxaux.resize(off);
    }
    // fill seeds
    {
        for(bidx = beg; bidx < end; bidx++) {
            bin = ref_kbmbinv(kbm->bins, bidx);
            if(bin->closed) continue;
            if(((u4i)bin->off + 1) * KBM_BIN_SIZE > kbm->reads->buffer[bin->ridx].rdlen)
                continue;
            off = kbm->reads->buffer[bin->ridx].rdoff + bin->off * KBM_BIN_SIZE;
            split_FIXP_kmers_kbm(kbm->rdseqs, off, KBM_BIN_SIZE, kbm->par->ksize,
                                 kbm->par->psize, kbm->par->kmer_mod, kmers);
            for(i = 0; i < 2; i++) {
                for(j = 0; j < kmers[i]->size; j++) {
                    f = ref_kmeroffv(kmers[i], j);
                    if(f->closed) continue;
                    u = get_kbmhash(kbm->hashs[0], f->kmer);
                    if(u == NULL || u->flt) continue;
                    x = ref_kbmkauxv(kbm->kauxs[0], offset_kbmhash(kbm->hashs[0], u));
                    kbm->bins->buffer[bidx].degree++;
                    if(x->cnt < u->tot) {
                        if(x->cnt && getval_bidx(kbm, x->off + x->cnt - 1) == bidx &&
                           kbm->vec_bidxaux[x->off + x->cnt - 1].dir == i) {
                            // repeated kmer within one bin
                        } else {
                            kbm->vec_bidxaux[x->off + x->cnt].bidx = bidx ;
                            kbm->vec_bidxaux[x->off + x->cnt].dir = i;
                            kbm->vec_bidxaux[x->off + x->cnt].koff = f->off >> 1;
                            x->cnt++;
                        }
                    }
                }
            }
        }
    }
    // sort seeds within a kmer
    {
        reset_iter_kbmhash(kbm->hashs[0]);
        while((u = ref_iter_kbmhash(kbm->hashs[0]))) {
            x = ref_kbmkauxv(kbm->kauxs[0], offset_kbmhash(kbm->hashs[0], u));
            // TODO DANGEROUS
            // clear_tmpbmerv(bms);
            // for(j = 0; j < x->cnt; j++) {
            //     push_tmpbmerv(bms, (kbm_tmp_bmer_t){getval_bidx(kbm, x->off + j),
            //                                         kbm->sauxs->buffer[x->off + j]});
            // }
            // dog_sort_array(bms->buffer, bms->size, kbm_tmp_bmer_t, num_cmpgt(a.bidx, b.bidx));
            // kbm->seeds->buffer[x->off + 0].bidx = bms->buffer[0].bidx & MAX_U4;
            // kbm->sauxs->buffer[x->off + 0].bidx = bms->buffer[0].bidx >> 32;
            // kbm->sauxs->buffer[x->off + 0] = bms->buffer[0].aux;
            // len = 1;
            // for(j = 1; j < x->cnt; j++) {
            //     if(bms->buffer[j].bidx < bms->buffer[j - 1].bidx) {
            //         continue;
            //     }
            //     kbm->seeds->buffer[x->off + len].bidx = bms->buffer[j].bidx & MAX_U4;
            //     kbm->sauxs->buffer[x->off + len].bidx = bms->buffer[0].bidx >> 32;
            //     kbm->sauxs->buffer[x->off + len] = bms->buffer[j].aux;
            //     len++;
            // }
            // x->cnt = len;
            auto beg_iter = kbm->vec_bidxaux.begin() + x->off;
            std::sort(beg_iter, beg_iter + x->cnt, [](auto& a, auto& b){return a.bidx < b.bidx;});
        }
    }
    free_kmeroffv(kmers[0]);
    free_kmeroffv(kmers[1]);
    // free_tmpbmerv(bms);
}

KBMDP *init_kbmdp() {
    KBMDP *dp;
    dp = (KBMDP *)malloc(sizeof(KBMDP));
    dp->kms = init_kbmdpev(1024);
    dp->km_len = 0;
    dp->cmask = init_bitvec(1024);
    dp->cms = init_kbmcmerv(64);
    dp->coffs = init_u4v(32);
    dp->rmask[0] = init_bitvec(256);
    dp->rmask[1] = init_bitvec(256);
    dp->cells[0] = init_kbmcellv(16);
    dp->cells[1] = init_kbmcellv(16);
    dp->bts = init_bit2vec(1024);
    dp->paths = init_kbmphash(13);
    dp->boff = 0;
    dp->last_bidx = 0;
    return dp;
}

void reset_kbmdp(KBMDP *dp, KBMAux *aux, u8i bidx) {
    //clear_kbmdpev(dp->kms);
    //dp->km_len = 0;
    clear_bitvec(dp->cmask);
    clear_kbmcmerv(dp->cms);
    clear_u4v(dp->coffs);
    recap_bitvec(dp->rmask[0], aux->qnbit);
    zeros_bitvec(dp->rmask[0]);
    recap_bitvec(dp->rmask[1], aux->qnbit);
    zeros_bitvec(dp->rmask[1]);
    recap_kbmcellv(dp->cells[0], aux->qnbin);
    recap_kbmcellv(dp->cells[1], aux->qnbin);
    clear_bit2vec(dp->bts);
    if(dp->paths->size > 1024 * 10) {
        free_kbmphash(dp->paths);
        dp->paths = init_kbmphash(1023);
    } else {
        clear_kbmphash(dp->paths);
    }
    dp->boff = bidx;
    dp->last_bidx = bidx;
}

void clear_kbmdp(KBMDP *dp, KBMAux *aux, u8i bidx) {
    reset_kbmdp(dp, aux, bidx);
    clear_kbmdpev(dp->kms);
    dp->km_len = 0;
}

void free_kbmdp(KBMDP *dp) {
    free_kbmdpev(dp->kms);
    free_bitvec(dp->cmask);
    free_kbmcmerv(dp->cms);
    free_u4v(dp->coffs);
    free_bitvec(dp->rmask[0]);
    free_bitvec(dp->rmask[1]);
    free_kbmcellv(dp->cells[0]);
    free_kbmcellv(dp->cells[1]);
    free_bit2vec(dp->bts);
    free_kbmphash(dp->paths);
    free(dp);
}

KBMAux *init_kbmaux(KBM *kbm) {
    KBMAux *aux;
    aux = (KBMAux *)malloc(sizeof(KBMAux));
    aux->kbm = kbm;
    aux->par = kbm->par;
    aux->qtag = NULL;
    aux->qseqs = NULL;
    aux->qsoff = 0;
    aux->qlen = 0;
    aux->slen = 0;
    aux->qidx = 0;
    aux->qnbin = 0;
    aux->qnbit = (aux->qnbin + 63) & 0xFFFFFFC0U;
    aux->bmin = 0;
    aux->bmax = MAX_U8;
    aux->koffs[0] = init_kmeroffv(32);
    aux->koffs[1] = init_kmeroffv(32);
    aux->refs = init_kbmrefv(64);
    aux->rank = init_u4v(64);
    aux->nheap = 1024;
    aux->heaps = (u4v **)malloc((aux->nheap + 1) * sizeof(u4v *));
    {
        u4i i;
        for(i = 0; i <= aux->nheap; i++) {
            aux->heaps[i] = init_u4v(8);
        }
    }
    aux->bmlen = 1;
    aux->bmoff = 0;
    aux->bmcnt = aux->kbm->bins->size;
    aux->binmap = NULL;
    aux->caches[0] = init_kbmdpev(64);
    aux->caches[1] = init_kbmdpev(64);
    aux->dps[0] = init_kbmdp();
    aux->dps[1] = init_kbmdp();
    aux->hits = init_kbmmapv(16);
    aux->cigars = init_bitsvec(1024, 3);
    aux->solids = NULL;
    aux->str = init_string(1024);
    return aux;
}

void free_kbmaux(KBMAux *aux) {
    free_kmeroffv(aux->koffs[0]);
    free_kmeroffv(aux->koffs[1]);
    free_kbmrefv(aux->refs);
    free_u4v(aux->rank);
    if(aux->heaps) {
        u4i i;
        for(i = 0; i <= aux->nheap; i++) {
            if(aux->heaps[i]) free_u4v(aux->heaps[i]);
        }
        free(aux->heaps);
    }
    if(aux->binmap) free(aux->binmap);
    free_kbmdpev(aux->caches[0]);
    free_kbmdpev(aux->caches[1]);
    free_kbmdp(aux->dps[0]);
    free_kbmdp(aux->dps[1]);
    free_kbmmapv(aux->hits);
    free_bitsvec(aux->cigars);
    free_string(aux->str);
    free(aux);
}

void print_exists_index_kbm(KBM *kbm, char *qtag, BaseBank *rdseqs, u8i seqoff,
                            u4i seqlen, kmeroffv *kmers[2], FILE *out) {
    KBMPar *par;
    kbm_kmer_t *u;
    kbm_kaux_t *x;
    kmer_off_t *f;
    u4i i, j;
    par = kbm->par;
    split_FIXP_kmers_kbm(rdseqs, seqoff, seqlen, par->ksize, par->psize, par->kmer_mod,
                         kmers);
    for(i = 0; i < 2; i++) {
        for(j = 0; j < kmers[i]->size; j++) {
            f = ref_kmeroffv(kmers[i], j);
            if(f->closed) continue;
            if(kbm->flags & 0x04) f->kidx = 0;
            u = get_kbmhash(kbm->hashs[f->kidx], f->kmer);
            if(u == NULL) {
                continue;
            }
            x = ref_kbmkauxv(kbm->kauxs[f->kidx], offset_kbmhash(kbm->hashs[f->kidx], u));
            fprintf(out, "%s\t%d\t%c\t0x%llx\t%u\t%u\t%c\n", qtag, f->off, "+-"[f->dir],
                    f -> kmer, x -> cnt, u -> tot, "YN"[u->flt]);
        }
    }
}

int _update_dp_path_kbm(KBMDP *dp, u8i end, kbm_cell_t *c) {
    int exists;
    kbm_path_t *p, P;
    P.beg = c->beg;
    p = prepare_kbmphash(dp->paths, P, &exists);
    if(exists) {
        if(p->score < c->score) {
            p->end = end;
            p->mat = c->mat;
            p->score = c->score;
        } else {
            return 0;
        }
    } else {
        p->beg = c->beg;
        p->end = end;
        p->mat = c->mat;
        p->score = c->score;
    }
    return 1;
}

void _dp_cal_spare_row_kbm(KBMAux *aux, int dir) {
    KBMDP *dp;
    kbmcellv *cells[2];
    BitVec *masks[2];
    kbm_cmer_t *m;
    kbm_cell_t D, H, V;
    u4i i, n, ni, rowoff;
    u8i bitoff, celoff;
    int flg, is_gap, score;
    dp = aux->dps[dir];
    flg = (dp->last_bidx - dp->boff) & 0x01;
    cells[0] = dp->cells[flg];
    cells[1] = dp->cells[!flg];
    masks[0] = dp->rmask[flg];
    masks[1] = dp->rmask[!flg];
    rowoff = dp->coffs->buffer[dp->last_bidx - dp->boff];
    bitoff = (dp->last_bidx - dp->boff) * aux->qnbit;
    celoff = (dp->last_bidx - dp->boff) * aux->qnbin;
    // update cells' mask to be calculated by new coming kbm_cmer_t
    for(i = 0; i < aux->qnbit; i += 64) {
        masks[1]->bits[i >> 6] =
            masks[0]->bits[i >> 6] | dp->cmask->bits[(i + bitoff) >> 6];
    }
    // dp core
    H = KBM_CELL_NULL;
    i = next_one_bitvec(masks[1], 0);
    D = (i && get_bitvec(masks[0], i - 1)) ? cells[0]->buffer[i - 1] : KBM_CELL_NULL;
    reg_zeros_bitvec(masks[0], 0, i);
    encap_bit2vec(dp->bts, aux->qnbin);
    n = 0;
    while(i < aux->qnbin) {
        if(get_bitvec(dp->cmask, i + bitoff)) {
            is_gap = 0;
            m = ref_kbmcmerv(dp->cms, rowoff + n);
            score = m->kmat + m->kcnt;
            n++;
        } else {
            is_gap = 1;
            m = (kbm_cmer_t *)&KBM_CMER_NULL;
            score = aux->par->pgap;
        }
        // horizontal
        {
            H.score += score;
            H.mat += m->kmat;
            //H.gap = 0;
        }
        if(H.var < 0) {
            H.score += -aux->par->pvar;    // score increases for abs(var) decreased
        } else {
            H.score += aux->par->pvar;
        }
        H.var++;
        H.bt = 1;    // horizontal backtrace, insertion for query sequence
        // diagonal
        {
            D.score += score;
            D.mat += m->kmat;
            //D.gap = 0;
        }
        D.var = 0;    // whether to reset var
        D.bt = 0;     // diagonal backtrace
        if(D.score >= H.score) {
            H = D;
        }
        // vertical
        V = D = get_bitvec(masks[0], i) ? cells[0]->buffer[i] : KBM_CELL_NULL;
        {
            V.score += score;
            V.mat += m->kmat;
            //V.gap = 0;
        }
        if(V.var > 0) {
            V.score += -aux->par->pvar;    // score increases for abs(var) decreased
        } else {
            V.score += aux->par->pvar;
        }
        V.var--;
        V.bt = 2;    // vertical backtrace
        if(V.score > H.score) {
            H = V;
        }
        if(is_gap) {
            H.gap++;
        } else {
            H.gap = 0;
        }
        set_bit2vec(dp->bts, dp->bts->size + i, H.bt);
        // init new path ID when there is no progenitor
        if(H.beg == 0) {
            H.beg = 1 + i + celoff;
        }
#if __DEBUG__
        if(KBM_LOG >= KBM_LOG_ALL) {
            fprintf(
                KBM_LOGF,
                "KBMLOG%d [x=%d, y=%llu, beg=%llu, score=%d, gap=%d, var=%d, bt=%d]\n",
                __LINE__, i, dp->last_bidx, (u8i)H.beg, H.score, H.gap, H.var, H.bt);
        }
#endif
        if(H.score > 0 && H.gap <= aux->par->max_bgap &&
           num_abs(H.var) <= aux->par->max_bvar) {
            // set cell, move next
            one_bitvec(masks[1], i);
            cells[1]->buffer[i] = H;
            if(is_gap == 0 && H.mat >= aux->par->min_mat)
                _update_dp_path_kbm(dp, 1 + i + celoff, &H);
            i++;
        } else if(D.score > 0) {
            // move next, may have diagonal hit
            zero_bitvec(masks[1], i);
            H = KBM_CELL_NULL;
            i++;
        } else {
            zero_bitvec(masks[1], i);
            ni = next_one_bitvec(masks[1], i + 1);
            if(i + 1 < ni) {
                reg_zeros_bitvec(masks[1], i + 1, ni);
                D = H = KBM_CELL_NULL;
            }
            i = ni;
        }
    }
    // prepare next row
    //NEXT_ROW:
    push_u4v(dp->coffs, dp->cms->size);
    zero_bitvec(dp->cmask, dp->cmask->n_bit);
    dp->cmask->n_bit += aux->qnbit;
    one_bitvec(dp->cmask, dp->cmask->n_bit);
    dp->bts->size += aux->qnbin;
    dp->last_bidx++;
}

int _backtrace_map_kbm(KBMAux *aux, int dir, kbm_path_t *p) {
    KBMDP *dp;
    kbm_map_t *hit;
    kbm_cmer_t *c;
    u8i cgoff, sidx;
    u4i i, mat, cnt, gap, cglen;
    int tmp, x, y, bt;
    dp = aux->dps[dir];
    hit = next_ref_kbmmapv(aux->hits);
    hit->qidx = aux->qidx;
    hit->qdir = dir;
    hit->tidx = aux->kbm->bins->buffer[dp->boff + p->beg / aux->qnbin].ridx;
    hit->tdir = 0;
    hit->qb = p->beg % aux->qnbin;
    hit->qe = p->end % aux->qnbin;
    hit->tb = p->beg / aux->qnbin;
    hit->te = p->end / aux->qnbin;
    hit->aln = num_min(hit->qe - hit->qb, hit->te - hit->tb);
    hit->aln++;
    if(hit->aln < aux->par->min_aln) {
        aux->hits->size--;
        return 0;
    }
    cgoff = aux->cigars->size;
    cglen = 0;
    cnt = 0;
    mat = 0;
    gap = 0;
    x = hit->qe;
    y = hit->te;
    while(x >= hit->qb && y >= hit->tb) {
        bt = get_bit2vec(dp->bts, x + y * aux->qnbin);
        if(get_bitvec(dp->cmask, x + y * aux->qnbit)) {
            c = ref_kbmcmerv(dp->cms, rank_bitvec(dp->cmask, x + y * aux->qnbit));
            cnt += c->kcnt;
            mat += c->kmat;
            push_bitsvec(aux->cigars, bt);
        } else {
            gap++;
            push_bitsvec(aux->cigars, 0x4 | bt);
        }
        switch(bt) {
            case 0:
                x--;
                y--;
                break;
            case 1: x--; break;
            default: y--; break;
        }
    }
    cglen = aux->cigars->size - cgoff;
    if(mat < (u4i)aux->par->min_mat ||
       mat < UInt(hit->aln * KBM_BSIZE * aux->par->min_sim) ||
       gap > (u4i)(hit->aln * aux->par->max_gap) || hit->aln < (int)aux->par->min_aln ||
       num_diff(hit->qe - hit->qb, hit->te - hit->tb) >
           (int)num_max(aux->par->aln_var * hit->aln, 1.0)) {
        aux->hits->size--;
        aux->cigars->size = cgoff;
        return 0;
    }
    if(aux->par->self_aln && aux->solids) {
        // Obsolete
        x = hit->qe;
        y = hit->te;
        while(x >= hit->qb && y >= hit->tb) {
            bt = get_bit2vec(dp->bts, x + y * aux->qnbin);
            if(get_bitvec(dp->cmask, x + y * aux->qnbit)) {
                c = ref_kbmcmerv(dp->cms, rank_bitvec(dp->cmask, x + y * aux->qnbit));
                for(i = 0; i < c->kcnt; i++) {
                    sidx = seed2solid_idx_kbm(aux->kbm, dp->kms->buffer + c->koff + i);
                    one_bitvec(aux->solids, sidx);    // Thread-unsafe, but no hurt
                }
            } else {
            }
            switch(bt) {
                case 0:
                    x--;
                    y--;
                    break;
                case 1: x--; break;
                default: y--; break;
            }
        }
    }
    hit->qe = (hit->qe + 1);
    hit->tb = (dp->boff + hit->tb - aux->kbm->reads->buffer[hit->tidx].binoff);
    hit->te = (dp->boff + hit->te - aux->kbm->reads->buffer[hit->tidx].binoff + 1);
    hit->mat = mat;
    hit->cnt = cnt;
    hit->gap = gap;
    hit->cgoff = cgoff;
    hit->cglen = cglen;
    //if(hit->qe > (int)aux->qlen) hit->qe = aux->qlen;
    //if(hit->te > (int)aux->kbm->reads->buffer[hit->tidx].rdlen) hit->te = aux->kbm->reads->buffer[hit->tidx].rdlen;
    if(dir) {
        tmp = aux->qnbin - hit->qb;
        hit->qb = aux->qnbin - hit->qe;
        hit->qe = tmp;
    }
#if __DEBUG__
    if(KBM_LOG) {
        fprintf(KBM_LOGF, "HIT\tQ[%d]\t%c\t%d\t%d", hit->qidx, "+-"[hit->qdir], hit -> qb,
                hit -> qe);
        fprintf(KBM_LOGF, "\tT[%d]\t%c\t%d\t%d", hit->tidx, "+-"[hit->tdir], hit -> tb,
                hit -> te);
        fprintf(KBM_LOGF, "\t%d\t%d\t%d\t%d\n", hit->mat, hit->aln, hit->cnt, hit->gap);
    }
#endif
    return 1;
}

int check_hit_cigar_kbm(kbm_map_t *hit, BitsVec *cigars) {
    u4i i, bt;
    int x, y;
    x = y = -1;
    i = hit->cglen;
    while(i) {
        bt = get_bitsvec(cigars, hit->cgoff + i - 1);
        bt = bt & 0x03;
        x += (0b0011 >> bt) & 0x01;
        y += (0b0101 >> bt) & 0x01;
        i--;
    }
    return !(x + 1 + hit->qb == hit->qe && y + 1 + hit->tb == hit->te);
}

void print_hit_kbm(KBM *kbm, char *qtag, u4i qlen, kbm_map_t *hit, BitsVec *cigars,
                   String *_str, FILE *out) {
    String *str;
    u8i coff;
    u4i clen, len, bt, _bt;
    if(hit->mat == 0) return;
    fprintf(out, "%s\t%c\t%d\t%d\t%d", qtag, "+-"[hit->qdir], qlen,
            hit -> qb *KBM_BIN_SIZE, hit -> qe *KBM_BIN_SIZE);
    fprintf(out, "\t%s\t%c\t%d\t%d\t%d", kbm->reads->buffer[hit->tidx].tag,
            "+-"[hit->tdir], kbm -> reads -> buffer[hit->tidx].rdlen,
            hit->tb *KBM_BIN_SIZE, hit->te *KBM_BIN_SIZE);
    fprintf(out, "\t%d\t%d\t%d\t%d\t", hit->mat, hit->aln * KBM_BIN_SIZE, hit->cnt,
            hit->gap);
    if(cigars) {
        str = _str ? _str : init_string(64);
        bt = len = 0;
        coff = hit->cgoff;
        clen = hit->cglen;
        clear_string(str);
        while(clen) {
            _bt = get_bitsvec(cigars, coff + clen - 1);
            if(_bt == bt) {
                len++;
            } else {
                if(len > 1) {
                    add_int_string(str, len);
                    add_char_string(str, "MID?mid?"[bt]);
                } else if(len == 1) {
                    add_char_string(str, "MID?mid?"[bt]);
                }
                bt = _bt;
                len = 1;
            }
            clen--;
        }
        if(len > 1) {
            add_int_string(str, len);
            add_char_string(str, "MID?mid?"[bt]);
        } else if(len == 1) {
            add_char_string(str, "MID?mid?"[bt]);
        }
        fputs(str->string, out);
        if(_str == NULL) free_string(str);
    } else {
        fputc('*', out);
    }
    fprintf(out, "\n");
}

void fprint_hit_kbm(KBMAux *aux, u4i hidx, FILE *out) {
    kbm_map_t *hit;
    hit = ref_kbmmapv(aux->hits, hidx);
    print_hit_kbm(aux->kbm, aux->qtag, aux->slen, hit, aux->cigars, aux->str, out);
}

void flip_hit_kbmaux(KBMAux *dst, KBMAux *src, u4i hidx) {
    kbm_map_t *h1, *h2;
    u4i t, i;
    h2 = next_ref_kbmmapv(dst->hits);
    h1 = ref_kbmmapv(src->hits, hidx);
    h2->qidx = h1->tidx;
    h2->qdir = h1->qdir;
    h2->tidx = h1->qidx;
    h2->tdir = h1->tdir;
    h2->cgoff = dst->cigars->size;
    h2->cglen = h1->cglen;
    h2->qb = h1->tb;
    h2->qe = h1->te;
    h2->tb = h1->qb;
    h2->te = h1->qe;
    h2->mat = h1->mat;
    h2->cnt = h1->cnt;
    h2->aln = h1->aln;
    h2->gap = h1->gap;
    if(h1->qdir) {
        for(i = 0; i < h1->cglen; i++) {
            t = get_bitsvec(src->cigars, h1->cgoff + h1->cglen - i - 1);
            if((t & 0x03)) {
                t = ((~(t & 0x03)) & 0x03) | (t & 0x4);
            }
            push_bitsvec(dst->cigars, t);
        }
    } else {
        for(i = 0; i < h1->cglen; i++) {
            t = get_bitsvec(src->cigars, h1->cgoff + i);
            if((t & 0x03)) {
                t = ((~(t & 0x03)) & 0x03) | (t & 0x4);
            }
            push_bitsvec(dst->cigars, t);
        }
    }
}

int _dp_path2map_kbm(KBMAux *aux, int dir) {
    KBMDP *dp;
    kbm_path_t *p;
    u4i ret;
    dp = aux->dps[dir];
    ret = 0;
    index_bitvec_core(dp->cmask, roundup_times(dp->cmask->n_bit, 64 * 8));
    reset_iter_kbmphash(dp->paths);
    while((p = ref_iter_kbmphash(dp->paths))) {
        p->beg--;    // 1-based to 0-based
        p->end--;
#if __DEBUG__
        if(KBM_LOG >= KBM_LOG_HIG) {
            fprintf(KBM_LOGF,
                    "KBMLOG%d\t%d\t%c\tkbm_path_t[%llu(%d:%d),%llu(%d:%d),%d]\n",
                    __LINE__, aux->qidx, "+-"[dir], (u8i)p -> beg,
                    (u4i)(p->beg % aux->qnbin), (u4i)((p->beg / aux->qnbin) + dp->boff),
                    (u8i)p -> end, (u4i)(p->end % aux->qnbin),
                    (u4i)((p->end / aux->qnbin) + dp->boff), p -> score);
        }
#endif
        if(_backtrace_map_kbm(aux, dir, p)) {
            ret++;
        }
    }
    //clear_kbmphash(dp->paths);
    return ret;
}

// KBM's tag2idx is wrongly loaded, need to be corrected
void rebuild_tag2idx_kbm(void *_kbm, size_t aux) {
    KBM *kbm;
    kbm_read_t *rd;
    u4i i;
    UNUSED(aux);
    kbm = (KBM *)_kbm;
    clear_cuhash(
        kbm->tag2idx);    // hash size is not changed, thus there won't have hash re-size
    for(i = 0; i < kbm->reads->size; i++) {
        rd = ref_kbmreadv(kbm->reads, i);
        if(rd->tag) put_cuhash(kbm->tag2idx, (cuhash_t){rd->tag, i});
    }
    kbm->flags |= 1LLU << 0;
}

int simple_chain_all_maps_kbm(kbm_map_t *srcs, u4i size, BitsVec *src_cigars,
                              kbm_map_t *dst, BitsVec *dst_cigars, float max_aln_var) {
    kbm_map_t *hit;
    u4i i, x, y, z, f;
    if(size < 2) return 0;
    sort_array(srcs, size, kbm_map_t, num_cmpgt(a.tb, b.tb));
    *dst = srcs[0];
    dst->cgoff = dst_cigars->size;
    append_bitsvec(dst_cigars, src_cigars, srcs[0].cgoff, srcs[0].cglen);
    for(i = 1; i < size; i++) {
        hit = srcs + i;
        if(dst->te > hit->tb) {
            goto FAILED;
        } else {
            y = hit->tb - dst->te;
            dst->te = hit->te;
        }
        if(dst->qdir) {
            if(hit->qe > dst->qb) {
                goto FAILED;
            } else {
                x = dst->qb - hit->qe;
                dst->qb = hit->qb;
            }
        } else {
            if(dst->qe > hit->qb) {
                goto FAILED;
            } else {
                x = hit->qb - dst->qe;
                dst->qe = hit->qe;
            }
        }
        dst->mat += hit->mat;
        dst->cnt += hit->cnt;
        dst->gap += hit->gap;
        dst->gap += num_max(x, y);
        z = num_min(x, y);
        f = 0x4 | 0;    // diagonal GAP
        pushs_bitsvec(dst_cigars, f, z);
        x -= z;
        y -= z;
        if(x > y) {
            z = x;
            f = 0x4 | 1;
        } else {
            z = y;
            f = 0x4 | 2;
        }
        pushs_bitsvec(dst_cigars, f, z);
        append_bitsvec(dst_cigars, src_cigars, hit->cgoff, hit->cglen);
    }
    dst->aln = num_min(dst->qe - dst->qb, dst->te - dst->tb);
    if(dst->aln * max_aln_var < num_diff(dst->qe - dst->qb, dst->te - dst->tb)) {
        goto FAILED;
    }
    dst->cglen = dst_cigars->size - dst->cgoff;
    return 1;
FAILED:
    dst_cigars->size = dst->cgoff;
    return 0;
}

size_t kbmreadv_deep_obj_desc_cnt(void *list, int idx) {
    if(idx == 0)
        return ((kbmreadv *)list)->size;
    else
        return 0;
}

size_t kbm_obj_desc_cnt(void *kbm, int idx) {
    UNUSED(kbm);
    if(idx == 8 || idx == 9)
        return KBM_N_HASH;
    else
        return 1;
}
const obj_desc_t kbmreadv_deep_obj_desc = {.tag = "kbmreadv_deep_obj_desc",
                                           .size = sizeof(kbmreadv),
                                           .n_child = 1,
                                           .mem_type = {1},
                                           .addr = {offsetof(kbmreadv, buffer)},
                                           .desc = {&kbm_read_t_obj_desc},
                                           .cnt = kbmreadv_deep_obj_desc_cnt,
                                           .post = NULL};
const obj_desc_t kbm_obj_desc = {
    .tag = "kbm_obj_desc",
    .size = sizeof(KBM),
    .n_child = 10,
    .mem_type = {1, 1, 1, 1, 1, 1, /*1, 1,*/ 2, 2},
    .addr = {offsetof(KBM, par), offsetof(KBM, rdseqs), offsetof(KBM, reads),
             offsetof(KBM, tag2idx), offsetof(KBM, bins), offsetof(KBM, binmarks),
             /*offsetof(KBM, seeds), offsetof(KBM, sauxs),*/ offsetof(KBM, hashs),
             offsetof(KBM, kauxs)},
    .desc = {&kbmpar_obj_desc, &basebank_obj_desc, &kbmreadv_deep_obj_desc,
             &cuhash_obj_desc, &kbmbinv_obj_desc, &bitvec_obj_desc, /*&kbmbmerv_obj_desc,
             &kbmbauxv_obj_desc,*/ &kbmhash_obj_desc, &kbmkauxv_obj_desc},
    .cnt = kbm_obj_desc_cnt,
    .post = rebuild_tag2idx_kbm};
// Please note that, kbm->tag2idx is not functional after mem_load, because we use cuhash_obj_desc instread of cuhash_deep_obj_desc

#include "kbm_opt.h"
