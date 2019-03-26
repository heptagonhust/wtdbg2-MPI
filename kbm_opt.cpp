#include "kbm_defines.h"
#include "debrain.h"
#include <thread>
#include <iostream>
#include <glog/logging.h>

volatile auto fake_init = []() {
    google::InitGoogleLogging("./wtdbg");
    return 0;
}();

const obj_desc_t kbm_read_t_obj_desc = {
    "kbm_read_t_obj_desc",       sizeof(kbm_read_t),     1,    {1},
    {offsetof(kbm_read_t, tag)}, {&OBJ_DESC_CHAR_ARRAY}, NULL, NULL};

void map_kbm(KBMAux *aux) {
#ifdef TEST_MODE
    if(aux->par->test_mode >= 4) return;
#endif
    KBM *kbm = aux->kbm;
    LOG(INFO) << "hptr" << aux->hptr << " bmlem" << aux->bmlen;
    for(; aux->hptr < aux->bmlen; aux->hptr++) {
        if(aux->hptr - aux->bmoff >= aux->nheap) {
            aux->bmoff += aux->nheap;
            for(int i = 0; i < aux->nheap; i++) {
                clear_u4v(aux->heaps[i]);
            }
            // for each aux->refs
            // easy gpu aux->refsy
            for(int i = 0; i < aux->refs->size; i++) {
                kbm_ref_t *ref = ref_kbmrefv(aux->refs, i);
                // if or bug
                while(ref->boff < ref->bend) {
                    int hidx = ref->bidx / aux->bmcnt;
                    int heap_id = hidx - aux->bmoff;
                    if(heap_id < aux->nheap) {
                        push_u4v(aux->heaps[heap_id], i);
                    }
                    break;
                }
            }
        }
        u4v *heap = aux->heaps[aux->hptr - aux->bmoff];
        if(heap->size == 0) {
            continue;
        }
        clear_kbmdpev(aux->caches[0]);
        clear_kbmdpev(aux->caches[1]);
        for(int i = 0; i < heap->size; i++) {
            int idx = heap->buffer[i];
            kbm_ref_t *ref = ref_kbmrefv(aux->refs, idx);
            while(1) {
                kbm_baux_t *saux = ref_kbmbauxv(kbm->sauxs, ref->boff);
                int pdir = (ref->dir ^ saux->dir);
                if(((aux->par->strand_mask >> pdir) & 0x01)) {
                    push_kbmdpev(aux->caches[pdir], (kbm_dpe_t){ref->poffs[pdir], idx,
                                                                ref->bidx, saux->koff});
                }
                ref->boff++;
                ref->bidx = getval_bidx(aux->kbm, ref->boff);
                if(ref->boff >= ref->bend) break;

                int hidx = ref->bidx / aux->bmcnt;
                if(hidx > aux->hptr) {
                    if(hidx - aux->bmoff < aux->nheap) {
                        push_u4v(aux->heaps[hidx - aux->bmoff], idx);
                    }
                    break;
                }
            }
        }
        NORMAL_AT_TEST(aux->par, 2) {
            if(aux->caches[0]->size * (aux->par->ksize + aux->par->psize) <
               UInt(aux->par->min_mat)) {
                aux->caches[0]->size = 0;
            } else {
                sort_array(aux->caches[0]->buffer, aux->caches[0]->size, kbm_dpe_t,
                           num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
            }
            if(aux->caches[1]->size * (aux->par->ksize + aux->par->psize) <
               UInt(aux->par->min_mat)) {
                aux->caches[1]->size = 0;
            } else {
                sort_array(aux->caches[1]->buffer, aux->caches[1]->size, kbm_dpe_t,
                           num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
            }
            //sort_array(aux->caches[0]->buffer, aux->caches[0]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
            //sort_array(aux->caches[1]->buffer, aux->caches[1]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
            // TODO: sort by bidx+koff is more reasonable, need to modify push_kmer_match_kbm too
        }
        NORMAL_AT_TEST(aux->par, 1) {
            for(int i = 0; i < 2; i++) {
                for(int j = 0; j < aux->caches[i]->size; j++) {
                    push_kmer_match_kbm(aux, i, aux->caches[i]->buffer + j);
                }
            }
            if(aux->hits->size >= aux->par->max_hit) return;
        }
    }
    if(aux->par->strand_mask & 0x01) push_kmer_match_kbm(aux, 0, NULL);
    if(aux->par->strand_mask & 0x02) push_kmer_match_kbm(aux, 1, NULL);
}
