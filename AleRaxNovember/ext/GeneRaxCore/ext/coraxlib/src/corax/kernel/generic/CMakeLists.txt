include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE core_clv.c
            core_clv.h
            core_derivatives.c
            core_derivatives.h
            core_likelihood.c
            core_likelihood.h
            core_pmatrix.c
            core_pmatrix.h
            fast_parsimony.c
            fast_parsimony.h
            weighted_parsimony.c
            weighted_parsimony.h
)
