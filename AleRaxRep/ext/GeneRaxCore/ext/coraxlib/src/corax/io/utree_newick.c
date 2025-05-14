#include <stdio.h>

#include "utree_io.h"

static char *newick_utree_recurse(const corax_unode_t *root,
                                  char *(*cb_serialize)(const corax_unode_t *),
                                  int level)
{
  char *newick;
  int   size_alloced = 0;
  assert(root != NULL);
  if (!root->next)
  {
    if (cb_serialize)
    {
      newick       = cb_serialize(root);
      size_alloced = (int)strlen(newick);
    }
    else
    {
      size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
    }
  }
  else
  {
    const corax_unode_t *start      = root->next;
    const corax_unode_t *snode      = start;
    char *               cur_newick = NULL;
    do {
      char *subtree =
          newick_utree_recurse(snode->back, cb_serialize, level + 1);
      if (subtree == NULL)
      {
        corax_set_error(CORAX_ERROR_MEM_ALLOC,
                        "Unable to allocate enough memory.");
        return NULL;
      }

      if (snode == start) { cur_newick = subtree; }
      else
      {
        assert(cur_newick);
        char *temp   = cur_newick;
        size_alloced = asprintf(&cur_newick, "%s,%s", temp, subtree);
        free(temp);
        free(subtree);
      }
      snode = snode->next;
    } while (snode != root);

    if (level > 0)
    {
      if (cb_serialize)
      {
        char *temp   = cb_serialize(root);
        size_alloced = asprintf(&newick, "(%s)%s", cur_newick, temp);
        free(temp);
      }
      else
      {
        size_alloced = asprintf(&newick,
                                "(%s)%s:%f",
                                cur_newick,
                                root->label ? root->label : "",
                                root->length);
      }
      free(cur_newick);
    }
    else
      newick = cur_newick;
  }

  if (size_alloced < 0)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "memory allocation during newick export failed");
    return NULL;
  }

  return newick;
}

char *utree_export_newick(const corax_unode_t *root,
                          int                  export_rooted,
                          double               root_brlen,
                          char *(*cb_serialize)(const corax_unode_t *))
{
  char *newick;
  char *subtree1;
  char *subtree2;
  int   size_alloced;

  if (!root) return NULL;

  if (!root->next) root = root->back;

  if (export_rooted)
  {
    assert(!cb_serialize);

    subtree1 = newick_utree_recurse(root->back, cb_serialize, 1);
    subtree2 = newick_utree_recurse(root, cb_serialize, 0);

    size_alloced = asprintf(&newick,
                            "(%s,(%s)%s:%f);",
                            subtree1,
                            subtree2,
                            root->label ? root->label : "",
                            root_brlen);
  }
  else
  {
    subtree1 = newick_utree_recurse(root->back, cb_serialize, 1);
    subtree2 = newick_utree_recurse(root, cb_serialize, 0);

    size_alloced = asprintf(&newick,
                            "(%s,%s)%s;",
                            subtree1,
                            subtree2,
                            root->label ? root->label : "");
  }

  free(subtree1);
  free(subtree2);

  if (size_alloced < 0)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "memory allocation during newick export failed");
    return NULL;
  }

  //  printf("newick: %s\n", newick);

  return (newick);
}

CORAX_EXPORT char *
corax_utree_export_newick(const corax_unode_t *root,
                          char *(*cb_serialize)(const corax_unode_t *))
{
  return utree_export_newick(root, 0, 0, cb_serialize);
}

CORAX_EXPORT char *corax_utree_export_newick_rooted(const corax_unode_t *root,
                                                    double root_brlen)
{
  return utree_export_newick(root, 1, root_brlen, NULL);
}
