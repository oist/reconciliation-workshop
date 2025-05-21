include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE common.c
            common.h
            parsimony.c
            parsimony.h
            partition.c
            partition.c
            repeats.c
            partition.h
            corax_assert.hpp
)
