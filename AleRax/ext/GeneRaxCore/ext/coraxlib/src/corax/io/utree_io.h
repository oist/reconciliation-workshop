/**
 * File containing the input and output operation for trees.
 */
#ifndef CORAX_IO_UTREE_IO_H_
#define CORAX_IO_UTREE_IO_H_

#include "corax/tree/utree.h"

/* attributes for corax_utree_show_ascii() */
#define CORAX_UTREE_SHOW_LABEL (1 << 0)
#define CORAX_UTREE_SHOW_BRANCH_LENGTH (1 << 1)
#define CORAX_UTREE_SHOW_CLV_INDEX (1 << 2)
#define CORAX_UTREE_SHOW_SCALER_INDEX (1 << 3)
#define CORAX_UTREE_SHOW_PMATRIX_INDEX (1 << 4)
#define CORAX_UTREE_SHOW_DATA (1 << 5)

/* structures for SVG visualization */

typedef struct corax_svg_attrib_s
{
  int    precision;
  long   width;
  long   font_size;
  long   tip_spacing;
  long   stroke_width;
  long   legend_show;
  long   legend_spacing;
  long   margin_left;
  long   margin_right;
  long   margin_bottom;
  long   margin_top;
  long   node_radius;
  double legend_ratio;
} corax_svg_attrib_t;


#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * This function will create a corax_utree_t from a newick _file_. If there is
   * an error, the function will return CORAX_ERROR (which happens to be 0). If a
   * file containing a rooted tree is passed to this function, an error will be
   * generated, and parsing will fail.
   *
   * @ingroup corax_utree_t
   */
  CORAX_EXPORT corax_utree_t *corax_utree_parse_newick(const char *filename);

  /**
   * This function will create a corax_utree_t from a newick _file_. If there is
   * an error, the function will return CORAX_ERROR (which happens to be 0). The
   * tree contained in the file is allowed to be rooted or unrooted.
   *
   * @ingroup corax_utree_t
   */
  CORAX_EXPORT corax_utree_t *
               corax_utree_parse_newick_rooted(const char *filename);

  /**
   * This function will create a corax_utree_t from a newick _file_. If there is
   * an error, the function will return CORAX_ERROR (which happens to be 0). If a
   * file containing a rooted tree is passed to this function, it will **unroot**
   * the tree, and return an unrooted tree.
   *
   * @ingroup corax_utree_t
   */
  CORAX_EXPORT corax_utree_t *
               corax_utree_parse_newick_unroot(const char *filename);

  /**
   * This function will create a corax_utree_t from a newick _string_. If there is
   * an error, the function will return CORAX_ERROR (which happens to be 0). If a
   * string containing a rooted tree is passed to this function, an error will be
   * generated, and parsing will fail.
   *
   * @ingroup corax_utree_t
   */
  CORAX_EXPORT corax_utree_t *corax_utree_parse_newick_string(const char *s);

  /**
   * This function will create a corax_utree_t from a newick _string_. If there is
   * an error, the function will return CORAX_ERROR (which happens to be 0). The
   * tree contained in the string is allowed to be rooted or unrooted.
   *
   * @ingroup corax_utree_t
   */
  CORAX_EXPORT corax_utree_t *
               corax_utree_parse_newick_string_rooted(const char *s);

  /**
   * @ingroup corax_utree_t
   * This function will create a corax_utree_t from a newick _string_. If there is
   * an error, the function will return CORAX_ERROR (which happens to be 0). If
   * the tree contained in the string is rooted, it will **unroot** the tree, and
   *
   * @ingroup corax_utree_t
   * return a rooted tree.
   */
  CORAX_EXPORT corax_utree_t *
               corax_utree_parse_newick_string_unroot(const char *s);

  /* functions in utree_newick.c */

  CORAX_EXPORT char *
  corax_utree_export_newick(const corax_unode_t *root,
                            char *(*cb_serialize)(const corax_unode_t *));

  CORAX_EXPORT char *corax_utree_export_newick_rooted(const corax_unode_t *root,
                                                      double root_brlen);

  /* functions in utree_ascii.c */

  CORAX_EXPORT void corax_utree_show_ascii(const corax_unode_t *tree,
                                           int                  options);

  /* functions in utree_svg.c */

  CORAX_EXPORT corax_svg_attrib_t *corax_svg_attrib_create(void);

  CORAX_EXPORT void corax_svg_attrib_destroy(corax_svg_attrib_t *attrib);

  CORAX_EXPORT int corax_utree_export_svg(corax_utree_t *           tree,
                                          corax_unode_t *           root,
                                          const corax_svg_attrib_t *attribs,
                                          const char *              filename);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_IO_UTREE_IO_H_ */
