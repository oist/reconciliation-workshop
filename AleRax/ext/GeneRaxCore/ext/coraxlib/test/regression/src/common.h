/*
    Copyright (C) 2015 Diego Darriba, Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
#ifndef COMMON_H_
#define COMMON_H_

#include "corax/corax.h"

extern const corax_state_t odd5_map[256];

/* parse attributes from the arguments */
unsigned int get_attributes(int argc, char **argv);

/* skip current test */
void         skip_test();

corax_partition_t *parse_msa(const char * filename,
                           unsigned int states,
                           unsigned int rate_cats,
                           unsigned int rate_matrices,
                           corax_utree_t *tree,
                           unsigned int attributes);

corax_partition_t *parse_msa_reduced(const char * filename,
                                   unsigned int states,
                                   unsigned int rate_cats,
                                   unsigned int rate_matrices,
                                   corax_utree_t *tree,
                                   unsigned int attributes,
                                   int max_sites);

/* callback function for traverse the utree */
int              cb_full_traversal(corax_unode_t *node);

/* displays a tree */
void show_tree (corax_unode_t * tree, int SHOW_ASCII_TREE);

/* print error and exit */
void  fatal(const char *format, ...) __attribute__((noreturn));
char *xstrdup(const char *s);
void *xmalloc(size_t size);

#endif /* COMMON_H_ */
