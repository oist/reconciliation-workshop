include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE binary_io_operations.c
            binary_io_operations.h
            binary.c
            binary.h
            fasta.c
            fasta.h
            maps.c
            maps.h
            newick.cpp
            newick.hpp
            output.c
            output.h
            phylip.c
            phylip.h
            utree_ascii.c
            utree_io.h
            utree_newick.c
            utree_svg.c
)
