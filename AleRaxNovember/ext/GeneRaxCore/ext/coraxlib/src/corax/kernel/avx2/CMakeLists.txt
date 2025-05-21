include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE core_clv_avx2.c
            core_clv_avx2.h
            core_derivatives_avx2.c
            core_derivatives_avx2.h
            core_likelihood_avx2.c
            core_likelihood_avx2.h
            core_pmatrix_avx2.c
            core_pmatrix_avx2.h
            fast_parsimony_avx2.c
            fast_parsimony_avx2.h
)
