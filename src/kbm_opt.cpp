#include "kbm_defines.h"
#include "debrain.h"
#include <thread>
#include <iostream>
#include "glog/logging.h"
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


void query_index_kbm(KBMAux *aux, char *qtag, u4i qidx, BaseBank *rdseqs, u8i seqoff,
                     u4i seqlen) {
    // kbm_kmer_t *u;
    // kbm_kaux_t *x;
    // kmer_off_t *f;
    u8i sidx, bmin, bmax;
    u4i hidx, next;
    u4i i, j, l, tot, mr, pdir;
    auto kbm = aux->kbm;
    auto par = aux->par;

    aux->qtag = qtag ? qtag : kbm->reads->buffer[qidx].tag;
    aux->qseqs = rdseqs;
    aux->qsoff = seqoff;
    aux->qlen = seqlen;
    aux->slen = (seqlen / KBM_BIN_SIZE) * KBM_BIN_SIZE;
    aux->qidx = qidx;
    aux->qnbin = aux->slen / KBM_BIN_SIZE;
    aux->qnbit = (aux->qnbin + 63) & 0xFFFFFFC0U;
    clear_kbmdp(aux->dps[0], aux, 0);
    clear_kbmdp(aux->dps[1], aux, 0);
    clear_kbmrefv(aux->refs);
    clear_u4v(aux->rank);
    clear_kbmdpev(aux->caches[0]);
    clear_kbmdpev(aux->caches[1]);
    clear_kbmmapv(aux->hits);
    clear_bitsvec(aux->cigars);
    RETURN_IF_TEST(par, 7);

    bmin = par->self_aln
               ? kbm->reads->buffer[qidx].binoff + kbm->reads->buffer[qidx].bincnt
               : 0;
    if(bmin < aux->bmin) bmin = aux->bmin;
    bmax = aux->bmax;
    if(par->self_aln == 2) {
        par->strand_mask = 2;    // 1: forward; 2: reverse; 3: both
    }
    split_FIXP_kmers_kbm(rdseqs, seqoff, aux->slen, par->ksize, par->psize, par->kmer_mod,
                         aux->koffs);
    RETURN_IF_TEST(par, 6);
    tot = 0;
    for(i = 0; i < 2; i++) {
        next = 0;
        for(j = 0; j < aux->koffs[i]->size; j++) {
            auto f = ref_kmeroffv(aux->koffs[i], j);
            if(f->closed) continue;
            if(kbm->flags & (1LLU << 2)) f->kidx = 0;
            auto u = get_kbmhash(kbm->hashs[f->kidx], f->kmer);
            if(u == NULL || u->flt || u->tot < par->kmin) {
                continue;
            }
            auto x = ref_kbmkauxv(kbm->kauxs[f->kidx], offset_kbmhash(kbm->hashs[f->kidx], u));
            auto ref = next_ref_kbmrefv(aux->refs);
            ref->mer = u;
            ref->aux = x;
            ref->kidx = f->kidx;
            ref->off = f->off;
            ref->dir = i;
            if(par->self_aln && aux->solids) {
                sidx = (kbm->reads->buffer[qidx].rdoff + ref->off) >> 1;
                ref->fine = get_bitvec(aux->solids, sidx);
            } else {
                ref->fine = 0;
            }
            ref->qbidx = ref->off / KBM_BIN_SIZE;
            ref->poffs[0] = ref->off;
            ref->poffs[1] = aux->slen - (ref->off + (aux->par->ksize + aux->par->psize));
            ref->boff = x->off;
            ref->bend = x->off + x->cnt;
            ref->bidx = getval_bidx(kbm, x->off);
            ref->closed = 0;
            {
                // Refine boundray
                if(par->self_aln == 2) {    // reverse complementary only
                    while(ref->boff < ref->bend &&
                          getval_bidx(kbm, ref->bend - 1) > bmin) {
                        ref->bend--;
                    }
                }
                while(ref->boff < ref->bend) {
                    if(ref->bidx < bmin) {
                        ref->boff++;
                        ref->bidx = getval_bidx(kbm, ref->boff);
                        continue;
                        //} else if(ref->b->bidx >= bmax){
                        //break;
                    }
                    break;
                }
                while(ref->bend > ref->boff) {
                    if(getval_bidx(kbm, ref->bend - 1) > bmax) {
                        ref->bend--;
                    } else {
                        break;
                    }
                }
                if(ref->boff >= ref->bend) {
                    ref->closed = 1;
                }
            }
            if(ref->closed) {
                aux->refs->size--;
                continue;
            }
            tot += x->cnt;
        }
    }
    
    if(par->self_aln && aux->solids) {
        // Obsolete
        sort_array(aux->refs->buffer, aux->refs->size, kbm_ref_t,
                   num_cmpgt(a.off, b.off));
        tot = 0;
        next = 0;
        for(i = 0; i < aux->refs->size; i++) {// I don't know what is  fine,solids and next
            auto ref = ref_kbmrefv(aux->refs, i);
            if(ref->closed) {
                continue;
            } else if(ref->fine) {
                tot += ref->bend - ref->boff;
                next = ref->off + (aux->par->ksize + aux->par->psize) / 2 + 1;
            } else if(ref->off >= next) {
                tot += ref->bend - ref->boff;
            } else {
                ref->boff = ref->bend;
                ref->closed = 1;
            }
        }
    } else if(aux->par->ksampling < KBM_BIN_SIZE && aux->refs->size) {
        sort_array(aux->refs->buffer, aux->refs->size, kbm_ref_t,
                   num_cmpgtx(a.qbidx, b.qbidx, b.bend - b.boff, a.bend - a.boff));
        tot = 0;
        for(i = j = 0; i < aux->refs->size; i++) {
            if(aux->refs->buffer[i].qbidx != aux->refs->buffer[j].qbidx) {
                if(aux->refs->buffer[j].qbidx) {    // skip the first and last bin
                    if((i - j) > aux->par->ksampling) {
                        l = j + aux->par->ksampling;
                        for(; j < l; j++) {
                            tot += aux->refs->buffer[j].bend - aux->refs->buffer[j].boff;
                        }
                        for(; j < i; j++) {
                            aux->refs->buffer[j].boff = aux->refs->buffer[j].bend;
                            aux->refs->buffer[j].bidx =
                                getval_bidx(kbm, aux->refs->buffer[j].boff);
                            aux->refs->buffer[j].closed = 1;
                        }
                    }
                }
                j = i;
            }
        }
        //sort_array(aux->refs->buffer, aux->refs->size, kbm_ref_t, num_cmpgt(a.off, b.off));
    }
    sort_array(aux->refs->buffer, aux->refs->size, kbm_ref_t, num_cmpgt(a.off, b.off));
    // estimate binmap
    aux->bmoff = 0;
    // ？？？？
    if(aux->refs->size) {
        mr = aux->par->min_mat / (aux->par->ksize + aux->par->psize);//200 / 
        if(mr < 512) mr = 512; //mr = 512
        aux->bmlen = tot / mr;
        if(aux->bmlen == 0) aux->bmlen = 1;
        aux->bmcnt = (aux->kbm->bins->size + aux->bmlen - 1) / aux->bmlen;
        if(aux->bmcnt < aux->qnbin * 50) {
            aux->bmcnt = aux->qnbin * 50;
            aux->bmlen = (aux->kbm->bins->size + aux->bmcnt - 1) / aux->bmcnt;
        }
        //fprintf("")
    } else {
        aux->bmlen = 1;
        aux->bmcnt = aux->kbm->bins->size;
    }
    if(0 && aux->bmlen > aux->nheap) {
        aux->bmlen = aux->nheap;
        aux->bmcnt = (aux->kbm->bins->size + aux->bmlen - 1) / aux->bmlen;
    }
    fprintf(stderr, " -- %s tot=%d bmlen=%d bmcnt=%d mr=%d in %s -- %s:%d --\n", aux->qtag, tot, aux->bmlen, aux->bmcnt, tot / aux->bmlen, __FUNCTION__, __FILE__, __LINE__); 
    for(i = 0; i < aux->nheap && i < aux->bmlen; i++) {
        clear_u4v(aux->heaps[i]);
    }
    //fprintf(stderr, " -- %s tot=%d avg=%d bmlen=%d bmcnt=%d mr=%d qnbin=%d in %s -- %s:%d --\n", aux->qtag, tot, tot / aux->bmlen, aux->bmlen, aux->bmcnt, mr, aux->qnbin, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
    RETURN_IF_TEST(par, 5);
    // init heaps
    for(i = 0; i < aux->refs->size; i++) {
        auto ref = ref_kbmrefv(aux->refs, i);
        while(ref->boff < ref->bend) {
            if(0) {
                pdir = (ref->dir ^ kbm->vec_bidxaux[ref->boff].dir);
                if(((aux->par->strand_mask >> pdir) & 0x01) == 0) {
                    ref->boff++;
                    ref->bidx = getval_bidx(aux->kbm, ref->boff);
                    continue;
                }
            }
            hidx = ref->bidx / aux->bmcnt;
            if(hidx - aux->bmoff < aux->nheap) {
                push_u4v(aux->heaps[hidx - aux->bmoff], i);
            }
            break;
        }
    }





    RETURN_IF_TEST(aux->par, 4);
    // auto kbm = aux->kbm;
    // auto par = aux->par;
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
                        if(aux->bmoff>0){
                            fprintf(stderr, " -- %s hidx=%d aux->bmoff=%d hidx - aux->bmoff=%d in %s -- %s:%d --\n", aux->qtag,hidx,aux->bmoff,hidx - aux->bmoff, __FUNCTION__, __FILE__, __LINE__);
                        }
                        push_u4v(aux->heaps[heap_id], i);
                    }
                }
            }
        }
        // if(aux->bmoff>0){
        //     fprintf(stderr, " -- %s aux->bmoff = %d in %s -- %s:%d --\n", aux->qtag,aux->bmoff, __FUNCTION__, __FILE__, __LINE__);
        // }
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
                    // fprintf(stderr, " -- %s hidx=%d aux->bmoff=%d hidx - aux->bmoff=%d in %s -- %s:%d --\n", aux->qtag,hidx,aux->bmoff,hidx - aux->bmoff, __FUNCTION__, __FILE__, __LINE__);
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

void map_kbm(KBMAux *aux) {
    // RETURN_IF_TEST(aux->par, 4);
    // // auto kbm = aux->kbm;
    // // auto par = aux->par;
    // LOG(INFO) << " bmlem" << aux->bmlen;

    // breakpoint();
    // for(int hptr = 0; hptr < aux->bmlen; hptr++) {
    //     if(hptr - aux->bmoff >= aux->nheap) {
    //         aux->bmoff += aux->nheap;
    //         for(int i = 0; i < aux->nheap; i++) {
    //             clear_u4v(aux->heaps[i]);
    //         }
    //         // for each aux->refs
    //         // easy gpu aux->refsy
    //         for(int i = 0; i < aux->refs->size; i++) {
    //             kbm_ref_t *kbm_ref = ref_kbmrefv(aux->refs, i);
    //             if(kbm_ref->boff < kbm_ref->bend) {
    //                 int hidx = kbm_ref->bidx / aux->bmcnt;
    //                 int heap_id = hidx - aux->bmoff;
    //                 if(heap_id < aux->nheap) {
    //                     if(aux->bmoff>0){
    //                         fprintf(stderr, " -- %s hidx=%d aux->bmoff=%d hidx - aux->bmoff=%d in %s -- %s:%d --\n", aux->qtag,hidx,aux->bmoff,hidx - aux->bmoff, __FUNCTION__, __FILE__, __LINE__);
    //                     }
    //                     push_u4v(aux->heaps[heap_id], i);
    //                 }
    //             }
    //         }
    //     }
    //     if(aux->bmoff>0){
    //         fprintf(stderr, " -- %s aux->bmoff = %d in %s -- %s:%d --\n", aux->qtag,aux->bmoff, __FUNCTION__, __FILE__, __LINE__);
    //     }
    //     u4v *heap = aux->heaps[hptr - aux->bmoff];

    //     if(heap->size == 0) {
    //         continue;
    //     }
    //     clear_kbmdpev(aux->caches[0]);
    //     clear_kbmdpev(aux->caches[1]);

    //     breakpoint();
    //     auto MASK = aux->par->strand_mask;
    //     for(int i = 0; i < heap->size; i++) {
    //         int idx = heap->buffer[i];
    //         kbm_ref_t *kbm_ref = ref_kbmrefv(aux->refs, idx);
    //         u8i max_bidx = aux->bmcnt * (hptr + 1);
    //         u8i boff = kbm_ref->boff; 
    //         #pragma unroll(8)
    //         for(;boff < kbm_ref->bend; ++boff) {
    //             // auto saux = ref_kbmbauxv(kbm->sauxs, boff);
    //             auto bidxaux = kbm->vec_bidxaux[boff];
    //             int pdir = (kbm_ref->dir ^ bidxaux.dir);
    //             if(MASK & (0x01 << pdir)) {
    //                 auto entry =
    //                     (kbm_dpe_t){kbm_ref->poffs[pdir], idx, kbm_ref->bidx, bidxaux.koff};
    //                 push_kbmdpev(aux->caches[pdir], entry);
    //             }
    //             kbm_ref->bidx = getval_bidx(kbm, boff + 1);
    //             if(kbm_ref->bidx >= max_bidx) break;
    //         }
    //         if(boff != kbm_ref->bend){
    //             ++boff;
    //             int hidx = kbm_ref->bidx / aux->bmcnt;
    //             if(hidx - aux->bmoff < aux->nheap) {
    //                 // fprintf(stderr, " -- %s hidx=%d aux->bmoff=%d hidx - aux->bmoff=%d in %s -- %s:%d --\n", aux->qtag,hidx,aux->bmoff,hidx - aux->bmoff, __FUNCTION__, __FILE__, __LINE__);
    //                 push_u4v(aux->heaps[hidx - aux->bmoff], idx);
    //             }
    //         }
    //         kbm_ref->boff = boff;
    //     }

    //     breakpoint();
    //     NORMAL_AT_TEST(par, 2) {
    //         for(int dir = 0; dir < 2; ++dir) {
    //             if(aux->caches[dir]->size * (par->ksize + par->psize) <
    //                UInt(par->min_mat)) {
    //                 aux->caches[dir]->size = 0;
    //             } else {
    //                 auto arr = aux->caches[dir]->buffer;
    //                 auto size = aux->caches[dir]->size;
    //                 // dog_sort_array(arr, size, kbm_dpe_t,
    //                 //            num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
    //                 std::sort(arr, arr + size, [](kbm_dpe_t &a, kbm_dpe_t &b) {
    //                     return a.bidx != b.bidx ? a.bidx < b.bidx : a.poff < b.poff;
    //                 });
    //             }
    //             // sort_array(aux->caches[dir]->buffer, aux->caches[dir]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
    //             // TODO: sort by bidx+koff is more reasonable, need to modify push_kmer_match_kbm too
    //         }
    //     }

    //     breakpoint();
    //     NORMAL_AT_TEST(par, 1) {
    //         for(int dir = 0; dir < 2; dir++) {
    //             for(int j = 0; j < aux->caches[dir]->size; j++) {
    //                 push_kmer_match_kbm(aux, dir, aux->caches[dir]->buffer + j);
    //             }
    //         }
    //         if(aux->hits->size >= par->max_hit) return;
    //     }
    // }
    // if(par->strand_mask & 0x01) push_kmer_match_kbm(aux, 0, NULL);
    // if(par->strand_mask & 0x02) push_kmer_match_kbm(aux, 1, NULL);
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
