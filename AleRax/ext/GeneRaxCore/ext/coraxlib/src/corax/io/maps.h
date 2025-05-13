#ifndef CORAX_IO_MAPS_H_
#define CORAX_IO_MAPS_H_

#include "corax/core/common.h"

CORAX_EXPORT extern const corax_state_t corax_map_bin[256];
CORAX_EXPORT extern const corax_state_t corax_map_nt[256];
CORAX_EXPORT extern const corax_state_t corax_map_aa[256];
CORAX_EXPORT extern const corax_state_t corax_map_gt10[256];
CORAX_EXPORT extern const unsigned int  corax_map_fasta[256];
CORAX_EXPORT extern const unsigned int  corax_map_phylip[256];
CORAX_EXPORT extern const unsigned int  corax_map_generic[256];

#endif /* CORAX_IO_MAPS_H_ */
