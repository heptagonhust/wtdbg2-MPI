#pragma once
#include "common.h"
#include "kbm.h"
EXTERN_C KBMPar *init_kbmpar();
EXTERN_C void free_kbmpar(KBMPar *par);
EXTERN_C KBM *init_kbm(KBMPar *par);
EXTERN_C void free_kbm(KBM *kbm);
EXTERN_C void reset_index_kbm(KBM *kbm);
EXTERN_C void clear_kbm(KBM *kbm);
EXTERN_C int cvt_kbm_read_length(u4i seqlen);
EXTERN_C void push_kbm(KBM *kbm, char *tag, int taglen, char *seq, u4i seqlen);
EXTERN_C void bitpush_kbm(KBM *kbm, char *tag, int taglen, u8i *seqs, u8i seqoff,
                          u4i seqlen);
EXTERN_C u8i filter_reads_kbm(KBM *kbm, u8i retain_size, int strategy);
EXTERN_C void ready_kbm(KBM *kbm);
EXTERN_C KBM *clone_seqs_kbm(KBM *src, KBMPar *par);
EXTERN_C void split_FIXP_kmers_kbm(BaseBank *rdseqs, u8i offset, u4i length, u1i ksize,
                                   u1i psize, u4i kmod, kmeroffv *rs[2]);
EXTERN_C u8i seed2solid_idx_kbm(KBM *kbm, kbm_dpe_t *p);
EXTERN_C u8i rdoff2solid_idx_kbm(KBM *kbm, u4i ridx, u4i roff);
EXTERN_C void index_kbm(KBM *kbm, u8i beg, u8i end, u4i ncpu, FILE *kmstat);
EXTERN_C void simple_index_kbm(KBM *kbm, u8i beg, u8i end);
EXTERN_C KBMDP *init_kbmdp();
EXTERN_C void reset_kbmdp(KBMDP *dp, KBMAux *aux, u8i bidx);
EXTERN_C void clear_kbmdp(KBMDP *dp, KBMAux *aux, u8i bidx);
EXTERN_C void free_kbmdp(KBMDP *dp);
EXTERN_C KBMAux *init_kbmaux(KBM *kbm);
EXTERN_C void free_kbmaux(KBMAux *aux);
EXTERN_C void query_index_kbm(KBMAux *aux, char *qtag, u4i qidx, BaseBank *rdseqs,
                              u8i seqoff, u4i seqlen);
EXTERN_C void print_exists_index_kbm(KBM *kbm, char *qtag, BaseBank *rdseqs, u8i seqoff,
                                     u4i seqlen, kmeroffv *kmers[2], FILE *out);
EXTERN_C int _update_dp_path_kbm(KBMDP *dp, u8i end, kbm_cell_t *c);
EXTERN_C void _dp_cal_spare_row_kbm(KBMAux *aux, int dir);
EXTERN_C int _backtrace_map_kbm(KBMAux *aux, int dir, kbm_path_t *p);
EXTERN_C int check_hit_cigar_kbm(kbm_map_t *hit, BitsVec *cigars);
EXTERN_C void print_hit_kbm(KBM *kbm, char *qtag, u4i qlen, kbm_map_t *hit,
                            BitsVec *cigars, String *_str, FILE *out);
EXTERN_C void fprint_hit_kbm(KBMAux *aux, u4i hidx, FILE *out);
EXTERN_C void flip_hit_kbmaux(KBMAux *dst, KBMAux *src, u4i hidx);
EXTERN_C int _dp_path2map_kbm(KBMAux *aux, int dir);
EXTERN_C void push_kmer_match_kbm(KBMAux *aux, int dir, kbm_dpe_t *p);
EXTERN_C void rebuild_tag2idx_kbm(void *_kbm, size_t aux);
EXTERN_C int simple_chain_all_maps_kbm(kbm_map_t *srcs, u4i size, BitsVec *src_cigars,
                                       kbm_map_t *dst, BitsVec *dst_cigars,
                                       float max_aln_var);
EXTERN_C size_t kbmreadv_deep_obj_desc_cnt(void *list, int idx);
EXTERN_C size_t kbm_obj_desc_cnt(void *kbm, int idx);