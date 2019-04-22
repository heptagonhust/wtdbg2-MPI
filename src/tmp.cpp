for(rid = qb; rid <= qe + ncpu; rid++) {
                if(rid < qe) {
                    // if(!KBM_LOG && ((rid - qb) % 2000) == 0) {
                    //     fprintf(KBM_LOGF, "\r%u|%llu", rid - qb, nhit);
                    //     fflush(KBM_LOGF);
                    // }
                    thread_wait_one(mdbg);
                } else {
                    thread_wait_next(mdbg);
                    pb = NULL;
                }
                if(mdbg->reg.closed == 0) {
                    KBMAux *aux = mdbg->aux;
                    if( g->corr_mode && mdbg->cc->cns->size) {
                        g->reads->buffer[mdbg->reg.rid].corr_bincnt =
                            mdbg->cc->cns->size / KBM_BIN_SIZE;
                    }
                    // if(alno) {
                    //     // beg_bufferedwriter(bw);
                    //     // if(g->corr_mode && mdbg->cc->cns->size) {
                    //     //     fprintf(bw->out, "#corrected\t%s\t%u\t",
                    //     //             mdbg->cc->tag->string, (u4i)mdbg->cc->cns->size);
                    //     //     println_fwdseq_basebank(mdbg->cc->cns, 0, mdbg->cc->cns->size,
                    //     //                             bw->out);
                    //     // }
                    //     for(i = 0; i < mdbg->aux->hits->size; i++) {
                    //         hit = ref_kbmmapv(mdbg->aux->hits, i);
                    //         fprint_hit_kbm(mdbg->aux, i, bw->out);
                    //     }
                    //     end_bufferedwriter(bw);
                    // }
                    for(i = 0; i < mdbg->aux->hits->size; i++) {
                        hit = ref_kbmmapv(mdbg->aux->hits, i);
                        if(hit->mat == 0) continue;
                        if(rdflags &&
                           g->kbm->reads->buffer[hit->tidx].bincnt <
                               g->kbm->reads->buffer[hit->qidx].bincnt &&
                           (hit->tb <= 1 &&
                            hit->te + 1 >=
                                (int)(g->kbm->reads->buffer[hit->tidx].bincnt)) &&
                           (hit->qb > 1 ||
                            hit->qe + 1 <
                                (int)(g->kbm->reads->buffer[hit->qidx].bincnt))) {
                            one_bitvec(rdflags, hit->tidx);
                        }
                    }
                    if(g->chainning_hits) {
                        chainning_hits_core(aux->hits, aux->cigars, g->uniq_hit,
                                            g->kbm->par->aln_var);
                    }
                    for(i = 0; i < aux->hits->size; i++) {
                        hit = ref_kbmmapv(aux->hits, i);
                        if(hit->mat == 0) continue;
                        //hit->qb  /= KBM_BIN_SIZE;
                        //hit->qe  /= KBM_BIN_SIZE;
                        //hit->tb  /= KBM_BIN_SIZE;
                        //hit->te  /= KBM_BIN_SIZE;
                        //hit->aln /= KBM_BIN_SIZE;
                        nhit++;
                        append_bitsvec(g->cigars, aux->cigars, hit->cgoff, hit->cglen);
                        hit->cgoff = g->cigars->size - hit->cglen;
                        if(raw) {
                            hit2rdregs_graph(
                                g, regs,
                                g->corr_mode ? mdbg->cc->cns->size / KBM_BIN_SIZE : 0,
                                hit, mdbg->aux->cigars, maps);
                        } else {
                            map2rdhits_graph(g, hit);
                        }
                    }
                    // if(KBM_LOG) {
                    //     fprintf(KBM_LOGF, "QUERY: %s\t+\t%d\t%d\n",
                    //             g->kbm->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg,
                    //             mdbg->reg.end);
                    //     for(i = 0; i < mdbg->aux->hits->size; i++) {
                    //         hit = ref_kbmmapv(mdbg->aux->hits, i);
                    //         fprintf(KBM_LOGF, "\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n",
                    //                 g->kbm->reads->buffer[hit->tidx].tag, "+-"[hit->qdir],
                    //                 g -> kbm -> reads -> buffer[hit->tidx].rdlen,
                    //                 hit->tb *KBM_BIN_SIZE, hit->te *KBM_BIN_SIZE,
                    //                 hit->aln *KBM_BIN_SIZE, hit->mat);
                    //     }
                    // }
                    mdbg->reg.closed = 1;
                }
                if(rid < qe && (rdflags == NULL || get_bitvec(rdflags, rid) == 0)) {
                    pb = ref_kbmreadv(g->kbm->reads, rid);
                    mdbg->reg = (reg_t){0, rid, 0, 0, pb->bincnt, 0, 0};
                    thread_wake(mdbg);
                }
            }
        }
        //if(!KBM_LOG) fprintf(KBM_LOGF, "\r%u reads|total hits %llu\n", qe - qb, nhit);
        if(reset_kbm) {
            reset_index_kbm(g->kbm);
        }
    }
    thread_beg_close(mdbg);
    free(mdbg->aux->par);
    free_kbmaux(mdbg->aux);
    if(g->corr_mode) {
        free_kbmblock((KBMBlock *)mdbg->cc->obj);
        free_ctgcns(mdbg->cc);
    }
    if(mdbg->raux) {
        free_kbm(mdbg->raux->kbm);
        free_kbmaux(mdbg->raux);
    }



