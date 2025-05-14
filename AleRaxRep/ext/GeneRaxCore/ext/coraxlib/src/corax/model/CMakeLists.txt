include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE aa_rates.c
            aa_rates.h
            eigen.c
            eigen.h
            evolmodel.cpp
            evolmodel.hpp
            gamma.c
            gamma.h
            invariant.c
            invariant.h
            modutil.c
            modutil.h
            modutil_aa.c
            modutil_dna.c
            modutil_gt.c
            modutil_mult.c
)
