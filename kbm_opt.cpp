#include "kbm_defines.h"
#include "debrain.h"
#include <thread>
#include <iostream>
#include <glog/logging.h>
#include <algorithm>
#include <atomic>

std::atomic_int gen_thread_idx(0);
volatile auto fake_init = []() {
    google::InitGoogleLogging("./wtdbg");
    return 0;
}();

void breakpoint() {
    thread_local auto thread_idx = gen_thread_idx++;
    if(thread_idx == 0) {
        ++fake_init;
    }
}

const obj_desc_t kbm_read_t_obj_desc = {
    "kbm_read_t_obj_desc",       sizeof(kbm_read_t),     1,    {1},
    {offsetof(kbm_read_t, tag)}, {&OBJ_DESC_CHAR_ARRAY}, NULL, NULL};

void map_kbm(KBMAux *aux) {
    RETURN_IF_TEST(aux->par, 4);
    auto kbm = aux->kbm;
    auto par = aux->par;
    LOG(INFO) << " bmlem" << aux->bmlen;

    breakpoint();
    for(int hptr = 0; hptr < aux->bmlen; hptr++) {
        if(hptr - aux->bmoff >= aux->nheap) {
            aux->bmoff += aux->nheap;
            for(int i = 0; i < aux->nheap; i++) {
                clear_u4v(aux->heaps[i]);
            }
            // for each aux->refs
            // easy gpu aux->refsy
            for(int i = 0; i < aux->refs->size; i++) {
                kbm_ref_t *kbm_ref = ref_kbmrefv(aux->refs, i);
                if(kbm_ref->boff < kbm_ref->bend) {
                    int hidx = kbm_ref->bidx / aux->bmcnt;
                    int heap_id = hidx - aux->bmoff;
                    if(heap_id < aux->nheap) {
                        push_u4v(aux->heaps[heap_id], i);
                    }
                }
            }
        }
        u4v *heap = aux->heaps[hptr - aux->bmoff];

        if(heap->size == 0) {
            continue;
        }
        clear_kbmdpev(aux->caches[0]);
        clear_kbmdpev(aux->caches[1]);

        breakpoint();
        auto MASK = aux->par->strand_mask;
        for(int i = 0; i < heap->size; i++) {
            int idx = heap->buffer[i];
            kbm_ref_t *kbm_ref = ref_kbmrefv(aux->refs, idx);
            u8i max_bidx = aux->bmcnt * (hptr + 1);
            u8i boff = kbm_ref->boff; 
            #pragma unroll(8)
            for(;boff < kbm_ref->bend; ++boff) {
                // auto saux = ref_kbmbauxv(kbm->sauxs, boff);
                auto bidxaux = kbm->vec_bidxaux[boff];
                int pdir = (kbm_ref->dir ^ bidxaux.dir);
                if(MASK & (0x01 << pdir)) {
                    auto entry =
                        (kbm_dpe_t){kbm_ref->poffs[pdir], idx, kbm_ref->bidx, bidxaux.koff};
                    push_kbmdpev(aux->caches[pdir], entry);
                }
                kbm_ref->bidx = getval_bidx(kbm, boff + 1);
                if(kbm_ref->bidx >= max_bidx) break;
            }
            if(boff != kbm_ref->bend){
                ++boff;
                int hidx = kbm_ref->bidx / aux->bmcnt;
                if(hidx - aux->bmoff < aux->nheap) {
                    push_u4v(aux->heaps[hidx - aux->bmoff], idx);
                }
            }
            kbm_ref->boff = boff;
        }

        breakpoint();
        NORMAL_AT_TEST(par, 2) {
            for(int dir = 0; dir < 2; ++dir) {
                if(aux->caches[dir]->size * (par->ksize + par->psize) <
                   UInt(par->min_mat)) {
                    aux->caches[dir]->size = 0;
                } else {
                    auto arr = aux->caches[dir]->buffer;
                    auto size = aux->caches[dir]->size;
                    // dog_sort_array(arr, size, kbm_dpe_t,
                    //            num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
                    std::sort(arr, arr + size, [](kbm_dpe_t &a, kbm_dpe_t &b) {
                        return a.bidx != b.bidx ? a.bidx < b.bidx : a.poff < b.poff;
                    });
                }
                // sort_array(aux->caches[dir]->buffer, aux->caches[dir]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
                // TODO: sort by bidx+koff is more reasonable, need to modify push_kmer_match_kbm too
            }
        }

        breakpoint();
        NORMAL_AT_TEST(par, 1) {
            for(int dir = 0; dir < 2; dir++) {
                for(int j = 0; j < aux->caches[dir]->size; j++) {
                    push_kmer_match_kbm(aux, dir, aux->caches[dir]->buffer + j);
                }
            }
            if(aux->hits->size >= par->max_hit) return;
        }
    }
    if(par->strand_mask & 0x01) push_kmer_match_kbm(aux, 0, NULL);
    if(par->strand_mask & 0x02) push_kmer_match_kbm(aux, 1, NULL);
}

void push_kmer_match_kbm(KBMAux *aux, int dir, kbm_dpe_t *p) {
    KBMDP *dp;
    kbm_cmer_t *c;
    kbm_dpe_t *e, E;
    u4i i, qb, bb, kmat, kcnt, blen;
    dp = aux->dps[dir];
    if(p == NULL) {
        if(dp->kms->size == 0) {
            dp->km_len = 0;
            return;
        }
        e = ref_kbmdpev(dp->kms, dp->kms->size - 1);
        if((int)dp->km_len < aux->par->min_mat ||
           e->bidx + 1 < dp->kms->buffer[0].bidx + aux->par->min_aln) {
            clear_kbmdpev(dp->kms);
            dp->km_len = aux->par->ksize + aux->par->psize;
            return;
        }
    } else {
        if(dp->kms->size == 0) {
            dp->km_len = (aux->par->ksize + aux->par->psize);
            push_kbmdpev(dp->kms, *p);
            return;
        }

        e = ref_kbmdpev(dp->kms, dp->kms->size - 1);
        if(e->bidx == p->bidx) {
            if(p->poff <= e->poff + aux->par->ksize + aux->par->psize) {
                dp->km_len += p->poff - e->poff;
            } else {
                dp->km_len += aux->par->ksize + aux->par->psize;
            }
            push_kbmdpev(dp->kms, *p);
            return;
        }

        if(e->bidx + aux->par->max_bgap + 1 < p->bidx ||
           dp->kms->buffer[0].bidx + aux->par->max_bcnt < p->bidx ||
           aux->kbm->bins->buffer[e->bidx].ridx != aux->kbm->bins->buffer[p->bidx].ridx) {
            if(Int(dp->km_len) < aux->par->min_mat ||
               e->bidx + 1 < dp->kms->buffer[0].bidx + aux->par->min_aln) {
                clear_kbmdpev(dp->kms);
                push_kbmdpev(dp->kms, *p);
                dp->km_len = aux->par->ksize + aux->par->psize;
                return;
            }
        } else {
            dp->km_len += aux->par->ksize + aux->par->psize;
            push_kbmdpev(dp->kms, *p);
            return;
        }
    }
    reset_kbmdp(dp, aux, dp->kms->buffer[0].bidx);

    // #ifdef TEST_MODE
    //     if(aux->par->test_mode >= 1) {
    //         clear_kbmdpev(dp->kms);
    //         if(p) {
    //             dp->km_len = aux->par->ksize + aux->par->psize;
    //             push_kbmdpev(dp->kms, *p);
    //         }
    //         return;
    //     }
    // #endif

    blen = dp->kms->buffer[dp->kms->size - 1].bidx + 1 - dp->kms->buffer[0].bidx + 2;
    {
        push_u4v(dp->coffs, 0);
        clear_bitvec(dp->cmask);
        encap_bitvec(dp->cmask, aux->qnbit * blen);
        reg_zeros_bitvec(dp->cmask, 0, aux->qnbit * blen);
        dp->cmask->n_bit = aux->qnbit;
        one_bitvec(dp->cmask, dp->cmask->n_bit);
    }
    qb = dp->kms->buffer[0].poff;
    bb = dp->kms->buffer[0].bidx;
    kcnt = 0;
    kmat = aux->par->ksize + aux->par->psize;
    E.bidx = dp->kms->buffer[dp->kms->size - 1].bidx + 1;
    E.poff = 0;
    for(i = 0; i <= dp->kms->size; i++) {
        e = (i < dp->kms->size) ? ref_kbmdpev(dp->kms, i) : &E;
        if(e->bidx == bb && e->poff / KBM_BIN_SIZE == qb / KBM_BIN_SIZE) {
            kcnt++;
            if(qb + aux->par->ksize + aux->par->psize >= e->poff) {
                kmat += e->poff - qb;
            } else {
                kmat += aux->par->ksize + aux->par->psize;
            }
        } else {
            one_bitvec(dp->cmask, (bb - dp->boff) * aux->qnbit + qb / KBM_BIN_SIZE);
            c = next_ref_kbmcmerv(dp->cms);
            c->koff = i - kcnt;
            c->kcnt = kcnt;
            c->kmat = kmat;
            c->boff = bb - dp->boff;
            while(bb < e->bidx) {
                _dp_cal_spare_row_kbm(aux, dir);
                bb++;
            }
            kcnt = 1;
            kmat = aux->par->ksize + aux->par->psize;
        }
        qb = e->poff;
    }
    // flush last row
    _dp_cal_spare_row_kbm(aux, dir);
    //collecting maps
    _dp_path2map_kbm(aux, dir);
    clear_kbmdpev(dp->kms);
    if(p) push_kbmdpev(dp->kms, *p);
    dp->km_len = aux->par->ksize + aux->par->psize;
}