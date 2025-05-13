#pragma once


#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <memory>

typedef struct corax_utree_s corax_utree_t;
typedef struct corax_unode_s corax_unode_t;
typedef struct corax_rtree_s corax_rtree_t;
typedef struct corax_rnode_s corax_rnode_t;
typedef unsigned long long corax_state_t;

char * corax_rtree_export_newick(const corax_rnode_t * root,
                                   char * (*cb_serialize)(const corax_rnode_t *));

struct PLLSequence {
  PLLSequence(char *label_, char *seq_, unsigned int len_):
    label(label_),
    seq(seq_),
    len(len_) {}
  char *label;
  char *seq;
  unsigned int len;
  ~PLLSequence() {
    free(label);
    free(seq);
  }
};
using PLLSequencePtr = std::unique_ptr<PLLSequence>;
using PLLSequencePtrs = std::vector<PLLSequencePtr>;
using SuperMatrix = std::unordered_map<std::string, std::string>;


class LibpllParsers {
public:
  LibpllParsers() = delete;

  static void parseMSA(const std::string &alignmentFilename, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights);

private:
  /**
   *  parse sequences and pattern weights from fasta file
   *  @param fasta_file Input file
   *  @param stateMap state std::map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parseFasta(const char *fasta_file, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights);

  /**
   *  parse sequences and pattern weights from phylip file
   *  @param phylip_file Input file
   *  @param stateMap state std::map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parsePhylip(const char *phylip_file, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights);
  


};

