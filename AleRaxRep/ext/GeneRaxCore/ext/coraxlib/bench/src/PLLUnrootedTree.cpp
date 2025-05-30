#include "PLLUnrootedTree.hpp"
#include "LibpllException.hpp"

#include <stack>
#include <functional>
#include <sstream>
#include <deque>
#include <algorithm>
#include <fstream>

void defaultUnodePrinter(corax_unode_t *node, 
    std::stringstream &ss)
{
  if (node->label) {
    ss << node->label;
  }
  ss << ":" << node->length;
}


static void destroyNodeData(void *)
{
}

void utreeDestroy(corax_utree_t *utree) {
  if(!utree)
    return;
  corax_utree_destroy(utree, destroyNodeData);
}



static corax_utree_t *readNewickFromStr(const std::string &str) 
{
  auto utree =  corax_utree_parse_newick_rooted(str.c_str());
  if (!utree) 
    throw LibpllException("Error while reading tree from std::string: ", str);
  return utree;
}


static corax_utree_t *readNewickFromFile(const std::string &str)
{
  try {
    auto utree =  corax_utree_parse_newick_rooted(str.c_str());
    if (!utree) 
      throw LibpllException("Error while reading tree from file: ", str);
    return utree;
  } catch (...) {
      throw LibpllException("Error while reading tree from file: ", str);
  }
}

static corax_utree_t *buildUtree(const std::string &str, bool isFile)
{
  if (isFile) {
    return readNewickFromFile(str);
  } else {
    return readNewickFromStr(str);
  }
}

PLLUnrootedTree::PLLUnrootedTree(const std::string &str, bool isFile):
  _tree(buildUtree(str, isFile), utreeDestroy)
{
}

std::unique_ptr<PLLUnrootedTree> PLLUnrootedTree::buildFromStrOrFile(const std::string &strOrFile)
{
  std::unique_ptr<PLLUnrootedTree> res;
  try {
    res = std::make_unique<PLLUnrootedTree>(strOrFile, true);
  } catch (...) {
    try {
      res = std::make_unique<PLLUnrootedTree>(strOrFile, false);
    } catch (...) {
    }
  }
  return res;
}



void PLLUnrootedTree::save(const std::string &fileName)
{
  std::ofstream os(fileName, std::ofstream::out);
  char *newick = corax_utree_export_newick_rooted(getRawPtr()->nodes[0], 0);
  os << newick;
  os.close();
  free(newick);
}

void PLLUnrootedTree::setMissingBranchLengths(double minBL)
{
  for (auto node: getLeaves()) {
    if (0.0 == node->length) {
      node->length = minBL;
    } 
  }
  for (unsigned int i = _tree->tip_count; i < _tree->tip_count + _tree->inner_count; ++i) {
    if (0.0 == _tree->nodes[i]->length)
      _tree->nodes[i]->length = minBL;
    if (0.0 == _tree->nodes[i]->next->length)
      _tree->nodes[i]->next->length = minBL;
    if (0.0 == _tree->nodes[i]->next->next->length)
      _tree->nodes[i]->next->next->length = minBL;
  }  
}
  
CArrayRange<corax_unode_t*> PLLUnrootedTree::getLeaves() const
{
  return CArrayRange<corax_unode_t*>(_tree->nodes, getLeavesNumber());
}

CArrayRange<corax_unode_t*> PLLUnrootedTree::getNodes() const
{
  return CArrayRange<corax_unode_t*>(_tree->nodes, getNodesNumber());
}

CArrayRange<corax_unode_t*> PLLUnrootedTree::getInnerNodes() const
{
  return CArrayRange<corax_unode_t*>(_tree->nodes + getLeavesNumber(), getInnerNodesNumber());
}


unsigned int PLLUnrootedTree::getNodesNumber() const
{
  return getLeavesNumber() + getInnerNodesNumber();
}

unsigned int PLLUnrootedTree::getDirectedNodesNumber() const
{
  return getLeavesNumber() + 3 * getInnerNodesNumber();
}

unsigned int PLLUnrootedTree::getLeavesNumber() const
{
  return _tree->tip_count;
}

unsigned int PLLUnrootedTree::getInnerNodesNumber() const
{
  return _tree->inner_count;
}
  
corax_unode_t *PLLUnrootedTree::getNode(unsigned int node_index) const
{
  return _tree->nodes[node_index];
}

corax_unode_t *PLLUnrootedTree::getAnyInnerNode() const
{
  return getNode(getLeavesNumber());
}
  
std::unordered_set<std::string> PLLUnrootedTree::getLeavesLabels()
{
  std::unordered_set<std::string> res;
  for (auto leaf: getLeaves()) {
    if (leaf->label) {
      res.insert(std::string(leaf->label));
    }
  }
  return res;
}


static bool isBranchIn(corax_unode_t *b, 
    const std::unordered_set<corax_unode_t *> &branches)
{
  return branches.find(b) != branches.end() 
    || branches.find(b->back) != branches.end();
}

std::unordered_set<corax_unode_t *> PLLUnrootedTree::getBranches() const
{
  std::unordered_set<corax_unode_t *> branches;
  for (auto node: getNodes()) {
    if (!isBranchIn(node, branches)) {
      branches.insert(node);
    }
    if (node->next) {
      node = node->next;
      if (!isBranchIn(node, branches)) {
        branches.insert(node);
      }
      node = node->next;
      if (!isBranchIn(node, branches)) {
        branches.insert(node);
      }
    }
  }
  return branches;
}

struct ToSort {
  corax_unode_t *node;
  std::string label;
  unsigned int distanceToLeaf;
 
  ToSort() {}
  ToSort(corax_unode_t *inode, 
      const std::string &ilabel,
      unsigned int idistanceToLeaf): node(inode),
  label(ilabel),
  distanceToLeaf(idistanceToLeaf) {}
  ToSort(const ToSort &s): node(s.node),
    label(s.label),
    distanceToLeaf(s.distanceToLeaf) {}
};

static void fillBranchesRec(corax_unode_t * node,
    std::vector<ToSort> &toSortVector,
    ToSort &toSort)
{
  if (!node->next) {
    toSort = ToSort(node, std::string(node->label), 0);
    toSortVector.push_back(toSort);
    return;
  }
  ToSort toSort1;
  fillBranchesRec(node->next->back, toSortVector, toSort1);
  ToSort toSort2;
  fillBranchesRec(node->next->next->back, toSortVector, toSort2);
  if (toSort1.label > toSort2.label) {
    toSort = ToSort(node, toSort1.label, toSort1.distanceToLeaf + 1);
  } else {
    toSort = ToSort(node, toSort2.label, toSort2.distanceToLeaf + 1);
  }
  toSortVector.push_back(toSort);
}

struct less_than_key
{
  inline bool operator() (const ToSort& t1, const ToSort& t2)
  {
    if (t1.label == t2.label) {
      return t1.distanceToLeaf < t2.distanceToLeaf;
    } 
    return t1.label < t2.label;
  }
};

std::vector<corax_unode_t *> PLLUnrootedTree::getBranchesDeterministic() const
{
  // find a deterministic root
  corax_unode_t *root = getNode(0);
  std::string rootLabel(root->label);
  for (auto leaf: getLeaves()) {
    std::string label(leaf->label);
    if (rootLabel > label) {
      rootLabel = label;
      root = leaf;
    }
  }
  std::vector<ToSort> toSortVector;
  ToSort toSort;
  fillBranchesRec(root->back, toSortVector, toSort);
  assert(toSortVector.size() == getLeavesNumber() * 2 - 3);
  std::sort(toSortVector.begin(), toSortVector.end(), less_than_key());
  std::vector<corax_unode_t *> res;
  for (const auto &toSort: toSortVector) {
    res.push_back(toSort.node);
  }
  return res;
}

static void fillPostOrder(corax_unode_t *node,
    std::vector<corax_unode_t*> &nodes,
    std::vector<char> &markedNodes)
{
  // we already traversed this node
  if (markedNodes[node->node_index]) {
    return;
  }
  // mark the node as traversed
  markedNodes[node->node_index] = true;
  // first process children
  if (node->next) {
    fillPostOrder(node->next->back, nodes, markedNodes);
    fillPostOrder(node->next->next->back, nodes, markedNodes);
  }
  nodes.push_back(node);
}
  
std::vector<corax_unode_t*> PLLUnrootedTree::getPostOrderNodesFrom(corax_unode_t *node) const
{
  std::vector<corax_unode_t*> nodes;
  std::vector<char> markedNodes(getDirectedNodesNumber(), false);
  fillPostOrder(node, nodes, markedNodes);
  return nodes;
}


std::vector<corax_unode_t*> PLLUnrootedTree::getPostOrderNodes(bool innerOnly) const
{
  std::vector<corax_unode_t*> nodes;
  std::vector<char> markedNodes(getDirectedNodesNumber(), false);
  if (innerOnly) {
    for (auto node: getLeaves()) {
      markedNodes[node->node_index] = true;
    }
  }
  // do the post order traversal from all possible virtual roots 
  for (auto node: getLeaves()) {
    fillPostOrder(node->back, nodes, markedNodes);
  }
  if (innerOnly) {
    assert(nodes.size() == getDirectedNodesNumber() - getLeavesNumber());
  } else {
    assert(nodes.size() == getDirectedNodesNumber());
  }
  return nodes;
}
  
std::vector<corax_unode_t*> PLLUnrootedTree::getReverseDepthNodes() const
{
  std::deque<corax_unode_t *> q;
  std::vector<bool> marked(getDirectedNodesNumber(), false);
  std::vector<corax_unode_t *> nodes;
  for (auto leaf: getLeaves()) {
    marked[leaf->node_index] = true;
    q.push_back(leaf);
  }
  while (!q.empty()) {
    auto node = q.front();
    q.pop_front();
    nodes.push_back(node);
    auto back = node->back;
    if (!back->next) {
      continue;
    }
    auto left = back->next;;
    auto right = back->next->next;
    if (!marked[left->node_index]) {
      q.push_back(left);
      marked[left->node_index] = true;
    }
    if (!marked[right->node_index]) {
      q.push_back(right);
      marked[right->node_index] = true;
    }
  }
  assert(nodes.size() == getDirectedNodesNumber());
  return nodes; 
}



static void getCladeRec(corax_unode_t *node, 
    std::unordered_set<unsigned int> &clade)
{
  if (node->next) {
    getCladeRec(node->next->back, clade);
    getCladeRec(node->next->next->back, clade);
  } else {
    clade.insert(node->node_index);
  }
}

std::unordered_set<unsigned int> 
  PLLUnrootedTree::getClade(corax_unode_t *node)
{
  std::unordered_set<unsigned int> clade;
  getCladeRec(node, clade);
  return clade;
}


// look for *pv (or one of his nexts) under the oriented node u.
// if it is found, *pv is updated with the found next, and 
// the function returns true (and false otherwise)
static bool orientAux(corax_unode_t *u, 
    corax_unode_t **pv,
    std::stack<corax_unode_t *> &path)
{
  path.push(u);
  auto v = *pv;
  if (v == u) {
    return true;
  }
  if (v->next) {
    if (v->next == u) {
      *pv = v->next;
      return true;
    } else if (v->next->next == u) {
      *pv = v->next->next;
      return true;
    }
  }
  if (!u->next) { 
    // end of recursion, we did not find *v
    path.pop();
    return false;
  } else {
    if (orientAux(u->next->back, pv, path) || 
      orientAux(u->next->next->back, pv, path)) {
      return true;
    } else {
      path.pop();
      return false;
    }
  }
}

static void stackToVector(std::stack<corax_unode_t *> s,
    std::vector<corax_unode_t *> &v)
{
  v.clear();
  v.resize(s.size());
  for (int i = s.size() - 1; i >= 0; --i) {
    v[i] = s.top();
    s.pop();
  }
}

void PLLUnrootedTree::orientTowardEachOther(corax_unode_t **pu,
    corax_unode_t **pv,
    std::vector<corax_unode_t *> &branchesPath)
{
  assert((*pu) != (*pv));
  auto *u = *pu;
  std::stack<corax_unode_t *> path;
  if (orientAux(u->back, pv, path)) {
    stackToVector(path, branchesPath);    
    return;
  }
  assert(path.size() == 0);
  if (u->next) {
    if (orientAux(u->next->back, pv, path)) {
      *pu = u->next;
      stackToVector(path, branchesPath);    
      return;
    } else if (orientAux(u->next->next->back, pv, path)) {
      *pu = u->next->next;
      stackToVector(path, branchesPath);    
      return;
    }
  }
  assert(false);
}




static void printAux(corax_unode_t *node,
    std::stringstream &ss,
    UnodePrinter f)
{
  if (node->next) {
    ss << "(";
    printAux(node->next->back, ss, f);
    ss << ",";
    printAux(node->next->next->back, ss, f);
    ss << ")";
  }
  f(node, ss);
}
  
std::string PLLUnrootedTree::getSubtreeString(corax_unode_t *subtree, UnodePrinter f)
{
  if (!subtree) {
    return "null";
  }
  std::stringstream ss;
  printAux(subtree, ss, f);
  return ss.str();
}

std::string PLLUnrootedTree::getNewickString(UnodePrinter f,
      corax_unode_t *root, 
      bool rooted)
{
  std::stringstream ss;
  if (!root) {
    root = getAnyInnerNode();
  }
  if (rooted) {
    ss << "(";
    printAux(root, ss, f);
    ss << ",";
    printAux(root->back, ss, f);
    ss << ");";
  } else {
    ss << "(";
    printAux(root->back, ss, f);
    ss << ",";
    printAux(root->next->back, ss, f);
    ss << ",";
    printAux(root->next->next->back, ss, f);
    ss << ");";
  }
  return ss.str();
}

std::unordered_set<std::string> PLLUnrootedTree::getLabels() const
{
  std::unordered_set<std::string> res;
  for (auto node:  getLeaves()) {
    if (node->label) {
      res.insert(node->label);
    }
  }
  return res;
}

static corax_unode_t *searchForSet(corax_unode_t *node, 
    std::unordered_set<std::string> &currentNodeSet,
    const std::unordered_set<std::string> &set)
{
  if (!node->next) {
    if (node->label) {
      currentNodeSet.insert(std::string(node->label));
    }
  } else {
    auto left = node->next->back;
    auto right = node->next->next->back;
    std::unordered_set<std::string> rightNodeSet;
    auto res1 = searchForSet(left, currentNodeSet, set);
    if (res1) {
      return res1;
    }
    auto res2 = searchForSet(right, rightNodeSet, set);
    if (res2) {
      return res2;
    }
    for (auto &elem: rightNodeSet) {
      currentNodeSet.insert(elem);
    }
  }
  return (currentNodeSet == set) ? node : nullptr;
}



static size_t leafHash(corax_unode_t *leaf) {
  assert(leaf);
  std::hash<std::string> hash_fn;
  return hash_fn(std::string(leaf->label));
}

static size_t getTreeHashRec(corax_unode_t *node, size_t i) {
  assert(node);
  if (i == 0) 
    i = 1;
  if (!node->next) {
    return leafHash(node);
  }
  auto hash1 = getTreeHashRec(node->next->back, i + 1);
  auto hash2 = getTreeHashRec(node->next->next->back, i + 1);
  std::hash<size_t> hash_fn;
  auto m = std::min(hash1, hash2);
  auto M = std::max(hash1, hash2);
  return hash_fn(m * i + M);

}

static corax_unode_t *findMinimumHashLeafRec(corax_unode_t * root, size_t &hashValue)
{
  assert(root);
  if (!root->next) {
    hashValue = leafHash(root);
    return root;
  }
  auto n1 = root->next->back;
  auto n2 = root->next->next->back;
  size_t hash1, hash2;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    hashValue = hash1;
    return min1;
  } else {
    hashValue = hash2;
    return min2;
  }
}

static corax_unode_t *findMinimumHashLeaf(corax_unode_t * root) 
{
  assert(root);
  auto n1 = root;
  auto n2 = root->back;
  size_t hash1 = 0;
  size_t hash2 = 0;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    return min1;
  } else {
    return min2;
  }
}

size_t PLLUnrootedTree::getUnrootedTreeHash() const
{
  auto minHashLeaf = findMinimumHashLeaf(getAnyInnerNode());
  auto res = getTreeHashRec(minHashLeaf, 0) + getTreeHashRec(minHashLeaf->back, 0);
  return res;
}



static const char *orderChildren(corax_unode_t *node,
    std::vector<bool> &leftFirst)
{
  if (!node->next) {
    return node->label;
  }
  const char *label1 = orderChildren(node->next->back, leftFirst);
  const char *label2 = orderChildren(node->next->next->back, leftFirst);
  if (strcmp(label1, label2) < 0) {
    leftFirst[node->node_index] = true;
    return label1;
  } else {
    leftFirst[node->node_index] = false;
    return label2;
  }
}

static bool areIsomorphicAux(corax_unode_t *node1,
    corax_unode_t *node2,
    const std::vector<bool> &leftFirst1,
    const std::vector<bool> &leftFirst2)
{
  if (!node1->next || !node2->next) {
    // at least one is a leaf
    if (node1->next || node2->next) {
      // only one is a leaf 
      return false;
    }
    return strcmp(node1->label, node2->label) == 0;
  }
  // both are internal nodes
  auto l1 = node1->next->back;
  auto r1 = node1->next->next->back;
  auto l2 = node2->next->back;
  auto r2 = node2->next->next->back;
  if (!leftFirst1[node1->node_index]) {
    std::swap(l1, r1);
  }
  if (!leftFirst2[node2->node_index]) {
    std::swap(l2, r2);
  }
  return areIsomorphicAux(l1, l2, leftFirst1, leftFirst2) 
    && areIsomorphicAux(r1, r2, leftFirst1, leftFirst2);
}

bool PLLUnrootedTree::areIsomorphic(const PLLUnrootedTree &t1,
    const PLLUnrootedTree &t2)
{
  if (t1.getNodesNumber() != t2.getNodesNumber()) {
    return false;
  }
  auto startingLeaf1 = findMinimumHashLeaf(t1.getAnyInnerNode())->back;
  auto startingLeaf2 = findMinimumHashLeaf(t2.getAnyInnerNode())->back;
  std::vector<bool> leftFirst1(t1.getDirectedNodesNumber(), true);;
  std::vector<bool> leftFirst2(t2.getDirectedNodesNumber(), true);;
  orderChildren(startingLeaf1, leftFirst1);
  orderChildren(startingLeaf2, leftFirst2);
  return areIsomorphicAux(startingLeaf1, 
      startingLeaf2, 
      leftFirst1, 
      leftFirst2);
}

bool PLLUnrootedTree::isBinary() const
{
  for (auto node: getInnerNodes()) {
    assert(node->next);
    if (node->next->next->next != node) {
      return false;
    }
  }
  return true;
}
  
corax_unode_t *PLLUnrootedTree::findLeaf(const std::string &label)
{
  for (auto leaf: getLeaves()) {
    if (label == leaf->label) {
      return leaf;
    }
  }
  return nullptr;
}

void PLLUnrootedTree::ensureUniqueLabels() 
{
  auto labels = getLabels();
  unsigned int i = 0;
  std::string prefix("s");
  for (auto node: getPostOrderNodes()) {
    if (node->next) {
      std::string newLabel;
      if (node->label) {
        newLabel = std::string(node->label);
      }
      while (labels.find(newLabel) != labels.end() || newLabel.size() == 0) {
        newLabel = prefix + std::to_string(i++);
      }
      free(node->label);
      node->label = static_cast<char*>(malloc(sizeof(char) * (newLabel.size() + 1)));
      std::strcpy(node->label, newLabel.c_str());
      labels.insert(newLabel);
    }
  }
}
  

