include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE bfgs.c
            brent.c
            callback.c
            callback.h
            em.c
            newton.c
            nni_round.c
            opt_branches.c
            opt_branches.h
            opt_generic.c
            opt_model.c
            opt_model.h
            opt_treeinfo.c
            opt_treeinfo.h
            spr_round.c
)

add_subdirectory(lbfgsb)
