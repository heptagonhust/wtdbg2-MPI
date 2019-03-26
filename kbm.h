/*
 *
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include "list.h"
#include "hashset.h"
#include "dna.h"
#include "filereader.h"
#include "bitvec.h"
#include "bitsvec.h"
#include "bit2vec.h"
#include "thread.h"
#include "txtplot.h"

//#define __DEBUG__	1
#define TEST_MODE

#define KBM_BSIZE 256
#define KBM_BIN_SIZE KBM_BSIZE
#define KBM_MAX_BINCNT 0xFFFFFFFFFFLLU    // 40 bits, 1024 G
#define KBM_MAX_RDCNT 0x3FFFFFFF          // 30 bits, 1 G
#define KBM_MAX_RDBINCNT 0xFFFFFF         // 24 bits
// able to index reference sequences
#define KBM_MAX_RDLEN 0xFFFFFFFFU    // 32 bits, 4 G bp

#define KBM_MAX_KSIZE 23
#define KBM_MAX_KCNT 0xFFFF    // 16 bits, 65535

#define KBM_N_HASH 4096

#define KBM_KF_BITS 32
#define KBM_KF_SIZE (1LLU << KBM_KF_BITS)
#define KBM_KF_MASK (KBM_KF_SIZE - 1LLU)

#define KBM_LOGF stderr
#define KBM_LOGFNO STDERR_FILENO
static int KBM_LOG = 0;
#define KBM_LOG_LOW 1
#define KBM_LOG_MID 2
#define KBM_LOG_HIG 3
#define KBM_LOG_ALL 4

#define KBM_MAX_RDGRP1 0x7FFFFF
#define KBM_MAX_RDGRP2 0xFF

typedef struct {
    u8i rdoff : 40, bincnt : 24;
    u4i rdlen, binoff;
    char *tag;
} kbm_read_t;
define_list(kbmreadv, kbm_read_t);

extern const obj_desc_t kbm_read_t_obj_desc;
extern const obj_desc_t kbmreadv_deep_obj_desc;



#if 0
#define KBM_MAX_BIN_DEGREE 0x7FFU
#endif
#define KBM_MAX_BIN_DEGREE 0x1FFU
// each BIN takes KBM_BIN_SIZE bp in uncompressed reads
typedef struct {
#if 0
	u4i ridx:28, off:24, closed:1, degree:11; // off * KBM_BIN_SIZE is the real position
#endif
    u4i ridx : 30, off : 24, closed : 1,
        degree : 9;    // off * KBM_BIN_SIZE is the real position
} kbm_bin_t;
define_list(kbmbinv, kbm_bin_t);

typedef struct {
    u4i bidx;
} kbm_bmer_t;
define_list(kbmbmerv, kbm_bmer_t);

typedef struct {
    u1i bidx;                 // bidx = (kbm_baux_t->bidx << 32 | kbm_bmer_t->bidx)
    u1i dir : 1, koff : 7;    // koff is the real (offset? >> 1), here offset is +0 or +1
} kbm_baux_t;
define_list(kbmbauxv, kbm_baux_t);

#define getval_bidx(kbm, offset)                          \
    ((((u8i)((kbm)->sauxs->buffer[offset].bidx)) << 32) | \
     (kbm)->seeds->buffer[offset].bidx)

//#define kbm_kmer_smear(K) ((K) ^ ((K) >> 4) ^ ((K) >> 7) ^ ((K) >> 12))
#define kbm_kmer_smear(K) ((K) ^ ((K) >> 4) ^ ((K) >> 7))

typedef struct {
    u8i mer : 46, tot : 17, flt : 1;
} kbm_kmer_t;
define_list(kbmkmerv, kbm_kmer_t);
#define KBM_KMERCODE(E) ((E).mer)
#define KBM_KMEREQUALS(E1, E2) ((E1).mer == (E2).mer)
#define KBM_KEYEQUALS(K, E) ((K) == (E).mer)
define_hashtable(kbmhash, kbm_kmer_t, KBM_KMERCODE, KBM_KMEREQUALS, u8i, ITSELF,
                 KBM_KEYEQUALS, kbm_kmer_t *, ITSELF);

typedef struct kbm_kaux_t{
    u8i off : 40, cnt : 24;
} kbm_kaux_t;
define_list(kbmkauxv, kbm_kaux_t);

typedef struct kbm_ref_t{
    kbm_kmer_t *mer;
    kbm_kaux_t *aux;
    u4i kidx;
    u4i off : 24, dir : 1, pdir : 1, fine : 1, closed : 1, extra_bits1 : 4;
    u4i qbidx;
    u4i poffs[2];
    u8i bidx, boff, bend;
    //kbm_bmer_t *b, *end;
} kbm_ref_t;
define_list(kbmrefv, kbm_ref_t);
#if 0
#define heap_cmp_kbm_bmer(refs, a, b)                                       \
    num_cmpx(refs[a].b->bidx, refs[b].b->bidx, refs[a].poffs[refs[a].pdir], \
             refs[b].poffs[refs[b].pdir])
#endif
#define heap_cmp_kbm_bmer(refs, a, b)                                   \
    num_cmpx(refs[a]->bidx, refs[b]->bidx, refs[a].poffs[refs[a].pdir], \
             refs[b].poffs[refs[b].pdir])

typedef struct {
    u4i koff;
    u4i kcnt : 8, kmat : 9, boff : 15;    // offset from the start bin_idx
} kbm_cmer_t;
define_list(kbmcmerv, kbm_cmer_t);
static const kbm_cmer_t KBM_CMER_NULL = {0, 0, 0, 0};

typedef struct {
    u8i beg : 46, mat : 16, bt : 2;
    b2i var;
    u2i gap;
    b4i score;
} kbm_cell_t;
static const kbm_cell_t KBM_CELL_NULL = {0, 0, 0, 0, 0, 0};
define_list(kbmcellv, kbm_cell_t);

typedef struct {
    u8i beg, end;
    u4i mat;
    int score;
} kbm_path_t;
define_list(kbmpathv, kbm_path_t);
#define kbmpath_hashcode(E) E.beg
#define kbmpath_hashequals(E1, E2) (E1).beg == (E2).beg
define_hashset(kbmphash, kbm_path_t, kbmpath_hashcode, kbmpath_hashequals);

typedef struct {
    u4i qidx : 31, qdir : 1;
    u4i tidx : 31, tdir : 1;
    u8i cgoff : 40, cglen : 24;
    int qb, qe, tb, te;
    int mat, cnt, aln, gap;    // gap is counted in BINs
} kbm_map_t;
define_list(kbmmapv, kbm_map_t);

typedef struct {
    int rd_len_order;    // 0
    //int hk; // 0
    int use_kf;                             // 0
    int min_bin_degree;                     // 2
    u4i ksize, psize;                       // 0, 21
    u4i kmax, kmin, kmer_mod, ksampling;    // 1000, 1, 4.0 * KBM_N_HASH, KBM_BSIZE
    float ktop;                             // 0.05
    // runtime
    u4i strand_mask;    // 3. 1: forward; 2: reverse; 3: both
    int self_aln;    // 0. 0: map to all; 1: only map to greater read_id; 2: itself but reverse complementary
    int skip_contained;      // 1
    u2i max_bgap;            // 4
    u2i max_bvar;            // 4
    float max_gap;           // 0.6
    u8i max_bcnt;            // 0xFFFF
    int pgap, pvar;          // -7, -21
    u4i max_hit;             // 1000
    int min_aln, min_mat;    // 2048/256, 200
    float aln_var;           // 0.25
    float min_sim;           // kmer similarity: 0.05
    int test_mode;           // see codes
} KBMPar;

static const obj_desc_t kbmpar_obj_desc = {
    "kbmpar_obj_desc", sizeof(KBMPar), 0, {}, {}, {}, NULL, NULL};

typedef struct {
    u8i flags;    // 64 bits, 0: mem_load all, 1: mem_load rdseqs+reads; 2: Single Hash Mode;3-63: unused
    KBMPar *par;
    BaseBank *rdseqs;
    kbmreadv *reads;
    cuhash *tag2idx;
    kbmbinv *bins;
    BitVec *binmarks;
    //u8i      *kfs;
    kbmbmerv *seeds;
    kbmbauxv *sauxs;
    kbmhash *hashs[KBM_N_HASH];
    kbmkauxv *kauxs[KBM_N_HASH];
} KBM;

extern const obj_desc_t kbm_obj_desc;

typedef struct {
#if 0
	u4i poff, bidx;
	u4i refidx:26, koff:6;
#endif
    u4i poff;
    u4i refidx;
    u8i bidx : 40, koff : 24;
} kbm_dpe_t;
define_list(kbmdpev, kbm_dpe_t);

typedef struct {
    kbmdpev *kms;    // kmer offset in query and bidx
    u4i km_len;
    BitVec *cmask;    // bit for kbm_cmer_t
    kbmcmerv *cms;
    u4v *coffs;    // kbm_cmer_t offset for each bin
    BitVec *rmask[2];
    kbmcellv *cells[2];
    Bit2Vec *bts;       // back trace flag: 0: diagonal, 1: horizontal, 2: vertical
    kbmphash *paths;    // storing best unique paths by now
    u8i boff;
    u8i last_bidx;
} KBMDP;

typedef struct {
    u8i kmer;
    u4i off;
    u4i kidx : 30, dir : 1, closed : 1;
} kmer_off_t;
define_list(kmeroffv, kmer_off_t);

typedef struct KBMAux {
    KBM *kbm;
    KBMPar *par;    // can diff from kbm->par
    char *qtag;
    BaseBank *qseqs;
    u8i qsoff;
    u4i qlen, slen, qnbin, qnbit;
    u4i qidx;
    u8i bmin, bmax;
    kmeroffv *koffs[2];
    kbmrefv *refs;
    u4v *rank, **heaps;
    u4i nheap;
    u4i hptr;
    u2i *binmap;
    u4i bmlen, bmoff, bmcnt;
    kbmdpev *caches[2];
    KBMDP *dps[2];
    kbmmapv *hits;
    BitsVec *cigars;
    BitVec *solids;
    String *str;
} KBMAux;


#include "kbm_opt.h"
#include "kbm_defines.h"