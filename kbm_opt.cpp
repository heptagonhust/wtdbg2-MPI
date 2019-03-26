#include "kbm_defines.h"
const obj_desc_t kbm_read_t_obj_desc = {
    "kbm_read_t_obj_desc",       sizeof(kbm_read_t),     1,    {1},
    {offsetof(kbm_read_t, tag)}, {&OBJ_DESC_CHAR_ARRAY}, NULL, NULL};



void map_kbm(KBMAux *aux) {
#ifdef TEST_MODE
    if(aux->par->test_mode >= 4) return;
#endif
    KBM *kbm = aux->kbm;
    while(aux->hptr < aux->bmlen) {
        if(aux->hptr - aux->bmoff >= aux->nheap) {
            aux->bmoff += aux->nheap;
            for(int i = 0; i < aux->nheap; i++) {
                clear_u4v(aux->heaps[i]);
            }
            for(int i = 0; i < aux->refs->size; i++) {
                kbm_ref_t *ref = ref_kbmrefv(aux->refs, i);
                while(ref->boff < ref->bend) {
                    int hidx = ref->bidx / aux->bmcnt;
                    if(hidx - aux->bmoff < aux->nheap) {
                        push_u4v(aux->heaps[hidx - aux->bmoff], i);
                    }
                    break;
                }
            }
        }
        u4v* heap = aux->heaps[aux->hptr - aux->bmoff];
        if(heap->size) {
            clear_kbmdpev(aux->caches[0]);
            clear_kbmdpev(aux->caches[1]);
            for(int i = 0; i < heap->size; i++) {
                int idx = heap->buffer[i];
                kbm_ref_t* ref = ref_kbmrefv(aux->refs, idx);
                while(1) {
                    kbm_baux_t *saux = ref_kbmbauxv(kbm->sauxs, ref->boff);
                    int pdir = (ref->dir ^ saux->dir);
                    if(((aux->par->strand_mask >> pdir) & 0x01)) {
                        push_kbmdpev(
                            aux->caches[pdir],
                            (kbm_dpe_t){ref->poffs[pdir], idx, ref->bidx, saux->koff});
                    }
                    ref->boff++;
                    ref->bidx = getval_bidx(aux->kbm, ref->boff);
                    if(ref->boff >= ref->bend) break;
#if __DEBUG__
                    if(ref->bidx < getval_bidx(aux->kbm, ref->boff - 1)) {
                        fprintf(stderr, " -- something wrong in %s -- %s:%d --\n",
                                __FUNCTION__, __FILE__, __LINE__);
                        fflush(stderr);
                        abort();
                    }
#endif
                    int hidx = ref->bidx / aux->bmcnt;
                    if(hidx > aux->hptr) {
                        if(hidx - aux->bmoff < aux->nheap) {
                            push_u4v(aux->heaps[hidx - aux->bmoff], idx);
                        }
                        break;
                    }
                }
            }
            {
#ifdef TEST_MODE
                if(aux->par->test_mode <= 2) {
#endif
                    if(aux->caches[0]->size * (aux->par->ksize + aux->par->psize) <
                       UInt(aux->par->min_mat)) {
                        aux->caches[0]->size = 0;
                    } else {
                        sort_array(aux->caches[0]->buffer, aux->caches[0]->size,
                                   kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
                    }
                    if(aux->caches[1]->size * (aux->par->ksize + aux->par->psize) <
                       UInt(aux->par->min_mat)) {
                        aux->caches[1]->size = 0;
                    } else {
                        sort_array(aux->caches[1]->buffer, aux->caches[1]->size,
                                   kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
                    }
                    //sort_array(aux->caches[0]->buffer, aux->caches[0]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
                    //sort_array(aux->caches[1]->buffer, aux->caches[1]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
                    // TODO: sort by bidx+koff is more reasonable, need to modify push_kmer_match_kbm too
#ifdef TEST_MODE
                }
#endif
#ifdef TEST_MODE
                if(aux->par->test_mode <= 1) {
#endif
                    for(int i = 0; i < 2; i++) {
                        for(int j = 0; j < aux->caches[i]->size; j++) {
#if __DEBUG__
                            if(KBM_LOG >= KBM_LOG_ALL) {
                                //fprintf(KBM_LOGF, "KBMLOG%d\t%d\t%d\t%c\t%d\t%llu[%d,%d]\t%llu[%d,%d]\n", __LINE__, aux->qidx, ref->poffs[ref->pdir], "+-"[ref->pdir], aux->hptr, ref->bidx, aux->kbm->bins->buffer[ref->bidx].ridx, aux->kbm->bins->buffer[ref->bidx].off * KBM_BIN_SIZE, (u8i)e->bidx, aux->kbm->bins->buffer[e->bidx].ridx, e->poff);
                            }
#endif
                            push_kmer_match_kbm(aux, i, aux->caches[i]->buffer + j);
                        }
                    }
                    if(aux->hits->size >= aux->par->max_hit) return;
#ifdef TEST_MODE
                }
#endif
            }
        }
        aux->hptr++;
    }
    if(aux->par->strand_mask & 0x01) push_kmer_match_kbm(aux, 0, NULL);
    if(aux->par->strand_mask & 0x02) push_kmer_match_kbm(aux, 1, NULL);
}


