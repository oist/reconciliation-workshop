Structures
================================================================================

A structure which simply keeps track of the "operations" required to calculate the likelihood on a tree, during the
partial likelihood calculation. In general, these should not be produced by hand, but instead by
`corax_utree_create_operations`.

```
typedef struct corax_operation
{
  unsigned int parent_clv_index;
  int parent_scaler_index;
  unsigned int child1_clv_index;
  unsigned int child1_matrix_index;
  int child1_scaler_index;
  unsigned int child2_clv_index;
  unsigned int child2_matrix_index;
  int child2_scaler_index;
} corax_operation_t;
```
