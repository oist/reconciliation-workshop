#include "LibpllParsers.hpp"
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <stack>
#include <array>
#include "LibpllException.hpp"
#include <corax/corax.h>
#include <iostream>






void LibpllParsers::parseMSA(const std::string &alignmentFilename, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights)
{
  if (!std::ifstream(alignmentFilename.c_str()).good()) {
    throw LibpllException("Alignment file " + alignmentFilename + "does not exist");
  }
  try {
    parseFasta(alignmentFilename.c_str(),
        stateMap, sequences, weights);
  } catch (...) {
    parsePhylip(alignmentFilename.c_str(),
        stateMap, sequences,
        weights);
  }
}

void LibpllParsers::parseFasta(const char *fastaFile, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights)
{
  auto reader = corax_fasta_open(fastaFile, corax_map_fasta);
  if (!reader) {
    throw LibpllException("Cannot parse fasta file ", fastaFile);
  }
  char * head;
  long head_len;
  char *seq;
  long seq_len;
  long seqno;
  int length;
  while (corax_fasta_getnext(reader, &head, &head_len, &seq, &seq_len, &seqno)) {
    sequences.push_back(PLLSequencePtr(new PLLSequence(head, seq, static_cast<unsigned int>(seq_len))));
    length = static_cast<int>(seq_len);
  }
  unsigned int count = static_cast<unsigned int>(sequences.size());
  char** buffer = static_cast<char**>(malloc(static_cast<size_t>(count) * sizeof(char *)));
  assert(buffer);
  for (unsigned int i = 0; i < count; ++i) {
    buffer[i] = sequences[i]->seq;
  }
  weights = corax_compress_site_patterns(buffer, stateMap, static_cast<int>(count), &length);
  if (!weights) 
    throw LibpllException("Error while parsing fasta: cannot compress sites from ", fastaFile);
  for (unsigned int i = 0; i < count; ++i) {
    sequences[i]->len = static_cast<unsigned int>(length);
  }
  free(buffer);
  corax_fasta_close(reader);
}
  
void LibpllParsers::parsePhylip(const char *phylipFile, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights)
{
  assert(phylipFile);
  assert(stateMap);
  std::unique_ptr<corax_phylip_t, void (*)(corax_phylip_t*)> reader(corax_phylip_open(phylipFile, corax_map_phylip),
      corax_phylip_close);
  if (!reader) {
    throw LibpllException("Error while opening phylip file ", phylipFile);
  }
  corax_msa_t *msa = nullptr;
  // todobenoit check memory leaks when using the std::exception trick
  try {
    
    msa = corax_phylip_parse_interleaved(reader.get());
    if (!msa) {
      throw LibpllException("failed to parse ", phylipFile);
    }
  } catch (...) {
    std::unique_ptr<corax_phylip_t, void(*)(corax_phylip_t*)> 
      reader2(corax_phylip_open(phylipFile, corax_map_phylip), corax_phylip_close);
    msa = corax_phylip_parse_sequential(reader2.get());
    if (!msa) {
      throw LibpllException("failed to parse ", phylipFile);
    }
  }
  weights = corax_compress_site_patterns(msa->sequence, stateMap, msa->count, &msa->length);
  if (!weights) 
    throw LibpllException("Error while parsing fasta: cannot compress sites");
  for (auto i = 0; i < msa->count; ++i) {
    PLLSequencePtr seq(new PLLSequence(msa->label[i], msa->sequence[i], static_cast<unsigned int>(msa->length)));
    sequences.push_back(std::move(seq));
    // avoid freeing these buffers with corax_msa_destroy
    msa->label[i] = nullptr;
    msa->sequence[i] = nullptr;
  }
  corax_msa_destroy(msa);
}
  

