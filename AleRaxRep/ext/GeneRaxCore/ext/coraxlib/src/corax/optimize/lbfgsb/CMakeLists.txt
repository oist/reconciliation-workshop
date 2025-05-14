include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET} PRIVATE lbfgsb.c lbfgsb.h linesearch.c linpack.c miniCBLAS.c
                                subalgorithms.c
)
