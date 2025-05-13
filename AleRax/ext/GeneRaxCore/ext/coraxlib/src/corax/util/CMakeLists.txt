include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE compress.c
            hardware.c
            list.c
            msa.c
            random.c
            stepwise.c
            absdiff.c
            nni_parsimony.c
    PUBLIC stepwise.h
           random.h
           msa.h
           list.h
           hardware.h
           compress.h
           absdiff.h
           nni_parsimony.h
)
