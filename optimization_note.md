# wtdbg2
1. sort_array / psort_array: hand-crafted sort func, seq & parallel. 
    - psort_array use C++ unsupported inline function definition, downgraded to seq one.
    - sort_array, use std::sort instead
    - std::sort with parallel policy is on the way
2. `sauxs`/`seeds` both hold one part of bidx. Setter/Getter, always treat them as one. 
    - merge `kbm_bmer_t` & `kbm_baux_t` into `kbm_bidx_t`. 
    - remove `kbm_temp_bmer_t`, which is served as tmp sorting array. 
    - disable mem_(load/find/dump)_obj. reconstruction of serialization will be considered later. 
3. `push_kbmdpev(aux->caches[pdir], (kbm_dpe_t){ref->poffs[pdir], idx, ref->bidx, saux->koff})`
    - seems to be the simliar reason 
    - investigate it later

