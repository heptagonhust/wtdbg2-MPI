#pragma once
#include "common.h"
#include "kbm.h"
KBMPar *init_kbmpar();
void free_kbmpar(KBMPar *par);
KBM *init_kbm(KBMPar *par);
void transfer_kbm(KBM *kbm, KBMPar *par,KBMPar *rpar,u4i *corr_mode,int world_rank);
void free_kbm(KBM *kbm);
void reset_index_kbm(KBM *kbm);
void clear_kbm(KBM *kbm);
int cvt_kbm_read_length(u4i seqlen);
void push_kbm(KBM *kbm, char *tag, int taglen, char *seq, u4i seqlen);
void bitpush_kbm(KBM *kbm, char *tag, int taglen, u8i *seqs, u8i seqoff,
                 u4i seqlen);
u8i filter_reads_kbm(KBM *kbm, u8i retain_size, int strategy);
void ready_kbm(KBM *kbm);
KBM *clone_seqs_kbm(KBM *src, KBMPar *par);
u8i seed2solid_idx_kbm(KBM *kbm, kbm_dpe_t *p);
u8i rdoff2solid_idx_kbm(KBM *kbm, u4i ridx, u4i roff);
void index_kbm(KBM *kbm, u8i beg, u8i end, u4i ncpu, FILE *kmstat);
void simple_index_kbm(KBM *kbm, u8i beg, u8i end);
KBMDP *init_kbmdp();
void reset_kbmdp(KBMDP *dp, KBMAux *aux, u8i bidx);
void clear_kbmdp(KBMDP *dp, KBMAux *aux, u8i bidx);
void free_kbmdp(KBMDP *dp);
KBMAux *init_kbmaux(KBM *kbm);
void free_kbmaux(KBMAux *aux);
void print_exists_index_kbm(KBM *kbm, char *qtag, BaseBank *rdseqs, u8i seqoff,
                            u4i seqlen, kmeroffv *kmers[2], FILE *out);
int _update_dp_path_kbm(KBMDP *dp, u8i end, kbm_cell_t *c);
void _dp_cal_spare_row_kbm(KBMAux *aux, int dir);
int _backtrace_map_kbm(KBMAux *aux, int dir, kbm_path_t *p);
int check_hit_cigar_kbm(kbm_map_t *hit, BitsVec *cigars);
void print_hit_kbm(KBM *kbm, char *qtag, u4i qlen, kbm_map_t *hit,
                   BitsVec *cigars, String *_str, FILE *out);
void fprint_hit_kbm(KBMAux *aux, u4i hidx, FILE *out);
void flip_hit_kbmaux(KBMAux *dst, KBMAux *src, u4i hidx);
int _dp_path2map_kbm(KBMAux *aux, int dir);
void push_kmer_match_kbm(KBMAux *aux, int dir, kbm_dpe_t *p);
void rebuild_tag2idx_kbm(void *_kbm, size_t aux);
int simple_chain_all_maps_kbm(kbm_map_t *srcs, u4i size, BitsVec *src_cigars,
                              kbm_map_t *dst, BitsVec *dst_cigars,
                              float max_aln_var);
size_t kbmreadv_deep_obj_desc_cnt(void *list, int idx);
size_t kbm_obj_desc_cnt(void *kbm, int idx);