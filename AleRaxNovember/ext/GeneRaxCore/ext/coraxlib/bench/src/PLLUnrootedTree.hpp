#pragma once

#include <corax/corax.h>
#include <unordered_set>
#include <string>
#include <memory>
#include <vector>
#include <unordered_set>
#include <functional>
#include <cstring>

/**
 * A wrapper arround C-style array to allow range iteration:
 *  int *array = new[size];
 *  ...
 *  for (int elem: CArrayRange<int>(array, size)) {
 *    // do something
 *  }
 */
template <typename DataType>
class CArrayRange {
public:
 class iterator {
 public:
   iterator(DataType * ptr): _ptr(ptr){}
   iterator operator++() { ++_ptr; return *this; }
   bool operator!=(const iterator & other) const { return _ptr != other._ptr; }
   const DataType& operator*() const { return *_ptr; }
 private:
   DataType* _ptr;
 };
private:
 DataType *_val;
 size_t _len;
public:
 CArrayRange(DataType *ptr, size_t len): _val(ptr), _len(len) {}
 iterator begin() const { return iterator(_val); }
 iterator end() const { return iterator(_val + _len); }
 size_t size() const { return _len; }
};



using UnodePrinter = 
  std::function<void(corax_unode_t *, std::stringstream &)>;

void defaultUnodePrinter(corax_unode_t *node, 
    std::stringstream &ss);

/**
 *  C++ wrapper around the libpll corax_utree_t structure
 *  to represent an unrooted tree
 */
class PLLUnrootedTree {
public:
  /**
   * Construct from a string that is either a path 
   * to a newick file or a newick string
   */
  PLLUnrootedTree(const std::string &str, bool isFile = true);



  static std::unique_ptr<PLLUnrootedTree> buildFromStrOrFile(const std::string &strOrFile);


  /**
   * Forbid copy
   */
  PLLUnrootedTree(const PLLUnrootedTree &) = delete;
  PLLUnrootedTree & operator = (const PLLUnrootedTree &) = delete;
  PLLUnrootedTree(PLLUnrootedTree &&) = delete;
  PLLUnrootedTree & operator = (PLLUnrootedTree &&) = delete;

  bool operator ==(const PLLUnrootedTree &other) const
  {
    return areIsomorphic(*this, other);
  }
  
  /*
   * Tree dimension
   */
  // number of nodes in the node buffer (#leaves + #innner)
  unsigned int getNodesNumber() const;
  // #leaves + 3 * #inner
  unsigned int getDirectedNodesNumber() const;
  // #leaves
  unsigned int getLeavesNumber() const;
  // #inner
  unsigned int getInnerNodesNumber() const; 
  
  std::unordered_set<std::string> getLabels() const;

  /*
   * Node access
   */
  corax_unode_t *getAnyInnerNode() const;
  corax_unode_t *getNode(unsigned int node_index) const;

 
  
  /**
   * labels
   */
  std::unordered_set<std::string> getLeavesLabels();

  /*
   * Save the tree in newick format in filename
   */
  void save(const std::string &fileName); 
  std::string getNewickString(UnodePrinter f = defaultUnodePrinter,
      corax_unode_t *root = nullptr, 
      bool rooted = false);
  static std::string getSubtreeString(corax_unode_t *subtree, UnodePrinter f = defaultUnodePrinter);

  /**
   * Replace null branch lengths with minBL
   */
  void setMissingBranchLengths(double minBL = 0.1); 

  
  size_t getUnrootedTreeHash() const;

  /**
   *  Direct access to the libpll structure
   */
  corax_utree_t *getRawPtr() {return _tree.get();}
  
  CArrayRange<corax_unode_t*> getLeaves() const;

  /*
   *  C++11 range for accessing internal nodes. 
   *  Only includes one of the three pll elements per node
   *  in the tree
   */
  CArrayRange<corax_unode_t*> getInnerNodes() const;

  /*
   *  C++11 range for accessing nodes. 
   *  For internal nodes, Only includes one of the three pll elements per node
   *  in the tree
   */
  CArrayRange<corax_unode_t*> getNodes() const;

  /**
   *  Create a vector of all nodes (including all three pll
   *  internal elements per internal node), such that a node
   *  always comes after its (virtual) children.
   */
  std::vector<corax_unode_t*> getPostOrderNodes(bool innerOnly = false) const;
  std::vector<corax_unode_t*> getPostOrderNodesFrom(corax_unode_t *node) const;
  std::vector<corax_unode_t*> getReverseDepthNodes() const;
  
 
  /**
   *  Return a set of all branches (for each node, one and only 
   *  one of node and node->back will be inserted in the set)
   */
  std::unordered_set<corax_unode_t *> getBranches() const;
  std::vector<corax_unode_t *> getBranchesDeterministic() const;


  /**
   *  replace pu and pv by one of their three possible next pointers
   *  such that their back pointers point to each other
   *  In addition, fill branchesPath with the sequence of BRANCHES
   *  along the path (if node1->back == node2, we only add one of the two)
   */
  static void orientTowardEachOther(corax_unode_t **pu, 
      corax_unode_t **pv,
      std::vector<corax_unode_t *> &branchesPath);



  /**
   *  Return the set of leaves under the input directed node
   */
  static std::unordered_set<unsigned int> getClade(corax_unode_t *node);

  static bool areIsomorphic(const PLLUnrootedTree &t1,
    const PLLUnrootedTree &t2);

  bool isBinary() const;
  
  void ensureUniqueLabels();

  /**
   *  Return the leaf node that has the label label,
   *  or null pointer if such leaf does not exist
   */
  corax_unode_t *findLeaf(const std::string &labe);

  static corax_unode_t *getLeft(corax_unode_t *node) {return node->next->back;}
  static corax_unode_t *getRight(corax_unode_t *node) {return node->next->next->back;}
private:
  std::unique_ptr<corax_utree_t, void(*)(corax_utree_t*)> _tree;
};





