corr_mode

hit->qidx
g->reads->buffer
g->corr_bsize
g->corr_bstep
g->corr_max
g->genome_size  
g->corr_gcov


```cpp
if(g->corr_mode) {
        mbp = g->genome_size * g->corr_gcov;
        qb = qe = g->kbm->reads->size / 2;
        nbp = g->kbm->reads->buffer[qb].rdlen;
        while(nbp < mbp && qb && qe + 1 < g->kbm->reads->size) {
            qb--;
            qe++;
            nbp += g->kbm->reads->buffer[qb].rdlen;
            nbp += g->kbm->reads->buffer[qe].rdlen;
        }
        if(qe < g->kbm->reads->size) qe++;
        fprintf(KBM_LOGF,
                "[%s] turn correct-mode on, reads[%u ~ %u = %u] (%llu bp), "
                "genome-size=%llu, corr-gcov=%0.2f, corr-dep=[%d,%d,%0.2f]\n",
                date(), qb, qe, qe - qb, nbp, g->genome_size, g->corr_gcov, g->corr_min,
                g->corr_max, g->corr_cov);
        fflush(KBM_LOGF);
    }

    if(g->corr_mode) {
        KBMBlock *kb;
        POGPar par;
        kb = init_kbmblock(g->corr_bsize, g->corr_bstep);
        //mdbg->cc = init_ctgcns(kb, iter_kbmblock, info_kbmblock, 1, 1, 1, g->corr_max, 200, 100, 1, 96, 2, -5, -2, -4, -1, 16, 3, 0.5, g->corr_bsize - g->corr_bstep + KBM_BIN_SIZE);
        par = DEFAULT_POG_PAR;
        par.refmode = 1;
        mdbg->cc =
            init_ctgcns(kb, iter_kbmblock, info_kbmblock, 1, 1, g->corr_max, 200, 100, 1,
                        g->corr_bsize - g->corr_bstep + KBM_BIN_SIZE, &par);
    }


in = g->corr_mode ? 1 : g->num_index;
if(g->kbm->vec_bidxaux.size()) {
    reset_kbm = 0;
    if(in > 1) {
        fprintf(KBM_LOGF, " ** WARNNING: change number of kbm index to 1 **\n");
        fflush(KBM_LOGF);
        in = 1;
    }
} else {
    reset_kbm = 1;
}
//fix_node = 0;
nhit = 0;
for(ii = 0; ii < in; ii++) {
    ib = ie;
    ie = ib + ic;
    while(ie > 0 && ie < g->kbm->bins->size &&
            g->kbm->bins->buffer[ie - 1].ridx == g->kbm->bins->buffer[ie].ridx)
        ie++;
    if(g->corr_mode == 0) {
        qb = 0;
        qe = ie ? g->kbm->bins->buffer[ie - 1].ridx : 0;
    }
    nbp = ((u8i)(ie - ib)) * KBM_BSIZE;
    if(reset_kbm) {
        reset_index_kbm(g->kbm);
        fprintf(KBM_LOGF, "[%s] indexing bins[%llu,%llu] (%llu bp), %d threads\n",
                date(), ib, ie, nbp, ncpu);
        fflush(KBM_LOGF);
        kmlog = (in > 1) ? NULL : open_file_for_write(prefix, ".kmerdep", 1);
        index_kbm(g->kbm, ib, ie, ncpu, kmlog);
        if(kmlog) {
            fclose(kmlog);
            kmlog = NULL;
        }
        fprintf(KBM_LOGF, "[%s] Done\n", date());
        fflush(KBM_LOGF);
        if(in == 1 && dump_kbm) {
            FILE *dump;
            fprintf(KBM_LOGF, "[%s] dump kbm index to %s ...", date(), dump_kbm);
            fflush(KBM_LOGF);
            dump = open_file_for_write(dump_kbm, NULL, 1);
            mem_dump_obj_file(g->kbm, 1, &kbm_obj_desc, 1, 0, dump);
            fclose(dump);
            fprintf(KBM_LOGF, " Done\n");
            fflush(KBM_LOGF);
        }
        if(in == 1) {
            u4i *deps;
            u8i hidx;
            kmlog = open_file_for_write(prefix, ".binkmer", 1);
            deps = (u4i*)calloc(KBM_BIN_SIZE + 1, 4);
            for(hidx = 0; hidx < g->kbm->bins->size; hidx++) {
                deps[g->kbm->bins->buffer[hidx].degree]++;
            }
            for(hidx = 0; hidx < KBM_BIN_SIZE; hidx++) {
                fprintf(kmlog, "%u\n", deps[hidx]);
            }
            fclose(kmlog);
            free(deps);
            if(!g->minimal_output) {
                kbm_bin_t *bn;
                kmlog = open_file_for_write(prefix, ".closed_bins", 1);
                for(hidx = 0; hidx < g->kbm->bins->size; hidx++) {
                    bn = ref_kbmbinv(g->kbm->bins, hidx);
                    if(bn->closed) {
                        fprintf(kmlog, "%s_F_%d_%d\t%d\n",
                                g->kbm->reads->buffer[bn->ridx].tag,
                                bn->off * KBM_BIN_SIZE, KBM_BIN_SIZE, bn->degree);
                    }
                }
                fclose(kmlog);
            }
        }
    }
```