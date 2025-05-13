include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE core_clv_sse.c
            core_clv_sse.h
            core_derivatives_sse.c
            core_derivatives_sse.h
            core_likelihood_sse.c
            core_likelihood_sse.h
            core_pmatrix_sse.c
            core_pmatrix_sse.h
            fast_parsimony_sse.c
            fast_parsimony_sse.h
)
