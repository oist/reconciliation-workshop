include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE core_clv_avx.c
            core_clv_avx.h
            core_derivatives_avx.c
            core_derivatives_avx.h
            core_likelihood_avx.c
            core_likelihood_avx.h
            core_pmatrix_avx.c
            core_pmatrix_avx.h
            fast_parsimony_avx.c
            fast_parsimony_avx.h
)
