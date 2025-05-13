#pragma once 

/*
 * This file was copied from RAxML-NG
 */

#include <string>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <random>


enum class FileFormat
{
  autodetect = 0,
  fasta,
  phylip,
  iphylip,
  vcf,
  catg,
  binary
};

enum class DataType
{
  autodetect = 0,
  dna,
  protein,
  binary,
  multistate,
  diploid10
};

enum class ParamValue
{
  undefined = 0,
  equal = 1,
  user = 2,
  model = 3,
  empirical = 4,
  ML = 5
};

enum class AscBiasCorrection
{
  none = 0,
  lewis = CORAX_ATTRIB_AB_LEWIS,
  felsenstein = CORAX_ATTRIB_AB_FELSENSTEIN,
  stamatakis = CORAX_ATTRIB_AB_STAMATAKIS,
};



typedef std::vector<double> doubleVector;
typedef std::vector<int> intVector;
typedef std::vector<unsigned int> uintVector;
typedef std::vector<std::string> NameList;
typedef std::pair<size_t,std::string> IdNamePair;
typedef std::vector<IdNamePair> IdNameVector;
typedef std::unordered_map<size_t,std::string> IdNameMap;
typedef std::unordered_map<std::string,size_t> NameIdMap;
typedef std::unordered_map<std::string,std::string> NameMap;
typedef std::set<size_t> IDSet;
typedef std::vector<size_t> IDVector;

typedef unsigned int WeightType;
typedef std::vector<WeightType> WeightVector;
typedef std::vector<WeightVector> WeightVectorList;
typedef std::unordered_map<size_t, WeightVector> WeightVectorMap;


/*
 * workaround needed for using enum as std::map key
 * code from: http://stackoverflow.com/a/24847480
 * */
struct EnumClassHash
{
  template <typename T>
  std::size_t operator()(T t) const
  {
      return static_cast<std::size_t>(t);
  }
};

