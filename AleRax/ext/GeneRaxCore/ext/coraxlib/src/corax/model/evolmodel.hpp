#ifndef CORAX_MODEL_EVOLMODEL_HPP_
#define CORAX_MODEL_EVOLMODEL_HPP_

#include <memory>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <stdexcept>

#include "corax/corax_core.h"
#include "modutil.h"

namespace corax {
namespace model {

typedef unsigned int corax_weight_t;
typedef std::unordered_map<corax_state_t,std::string> StateNameMap;

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

enum class DataType
{
  autodetect = 0,
  dna,
  protein,
  binary,
  multistate,
  genotype10
};

/**
 * Model parameter estimation mode. Possible values:
 * - undefined = not specified / default, see EvolModel::set_param_mode_default()
 * - equal     = all values are equal (e.g. all substitution rates are 1.0)
 * - user      = fixed user-defined values
 * - model     = fixed values specified by the pre-defined substitution matrix (e.g. LG)
 * - empirical = counted from the alignment (e.g. state frequencies)
 * - ML        = estimated by maximum likelihood from the alignment
 */
enum class ParamMode
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

/**
 * Time-reversible substitution matrix (Q-matrix) defined by two components:
 * 1. Vector of substitution rates ( = upper triangle of the symmetric rate matrix)
 * 2. Vector of stationary base frequencies
 */
class SubstitutionMatrix
{
public:
  SubstitutionMatrix(const corax_subst_model_t& sm) :
    _states(sm.states), _name(sm.name)
  {
    if (sm.freqs)
      _base_freqs.assign(sm.freqs, sm.freqs + sm.states);
    if (sm.rates)
      _subst_rates.assign(sm.rates, sm.rates + sm.states*(sm.states-1)/2);
    if (sm.rate_sym)
      _rate_sym.assign(sm.rate_sym, sm.rate_sym + sm.states*(sm.states-1)/2);
    if (sm.freq_sym)
      _freq_sym.assign(sm.freq_sym, sm.freq_sym + sm.states);
  };

  // getters

  /**
   * @brief Number of states (=matrix dimension), e.g. 4 for DNA
   */
  unsigned int states() const;

  /**
   * @brief Matrix name, e.g. GTR or LG
   */
  std::string name() const;

  /**
   * @brief Return array of substitution rates. Dimension: num_rates()
   */
  const std::vector<double>& subst_rates() const { return _subst_rates; }

  /**
   * @brief Return array of stationary base frequency. Dimension: states()
   */
  const std::vector<double>& base_freqs() const { return _base_freqs; }

  /**
   * @brief Return symmetry vector for the substitution rates. e.g.
   * GTR = 0 1 2 3 4 5
   * HKY = 0 1 0 0 1 0
   * JC  = 0 0 0 0 0 0
   */
  const std::vector<int>& rate_sym() const { return _rate_sym; }

  /**
   * @brief Return symmetry vector for the base frequencies
   */
  const std::vector<int>& freq_sym() const { return _freq_sym; }

  /**
   * @brief Return number of substitution rates in the matrix
   */
  unsigned int num_rates() const  { return _states*(_states-1)/2; }

  /**
   *  @brief Return number of *unique* rates considering rate symmetries,
   *  e.g. 6 for GTR and 2 for HKY
   */
  unsigned int num_uniq_rates() const
  {
    if (_rate_sym.empty())
      return num_rates();
    else
      return *std::max_element(_rate_sym.cbegin(), _rate_sym.cend()) + 1u;
  }

  /**
   * @brief Return array of all *unique* substitution rates, see also: num_uniq_rates()
   */
  std::vector<double> uniq_subst_rates() const
  {
    if (!_rate_sym.empty())
    {
      std::vector<double> uniq_rates(num_uniq_rates());
      for (size_t i = 0; i < _subst_rates.size(); ++i)
        uniq_rates[_rate_sym[i]] = _subst_rates[i];
      return uniq_rates;
    }
    else
      return _subst_rates;
  };


  // setters - see documentation for the respective getters above
  void base_freqs(const std::vector<double>& v)
  {
//    std::cout << "expected: " << _states << ", got: " << v.size() << std::endl;
    if (v.size() != _states)
      throw std::invalid_argument("Invalid size of base_freqs vector!");

    _base_freqs = v;
  };

  void subst_rates(const std::vector<double>& v)
  {
    if (v.size() != num_rates())
      throw std::invalid_argument("Invalid size of subst_rates vector!");

    _subst_rates = v;
  };

  void uniq_subst_rates(const std::vector<double>& v)
  {
    if (!_rate_sym.empty())
    {
      if (v.size() != num_uniq_rates())
        throw std::invalid_argument("Invalid size of subst_rates vector!");

      _subst_rates.resize(num_rates());
      for (size_t i = 0; i < _subst_rates.size(); ++i)
        _subst_rates[i] = v[_rate_sym[i]];
    }
    else
      subst_rates(v);
  };

private:
  unsigned int _states;
  std::string _name;
  std::vector<double> _base_freqs;
  std::vector<double> _subst_rates;
  std::vector<int> _rate_sym;
  std::vector<int> _freq_sym;
};

/**
 * Evolutionary model representation inspired by raxml-ng. Some features:
 * - parsing from/printing to string
 * - parsing substitution rates/freqs from PAML file
 * - parameter sync with corax_partition_t structs
 * - built-in DNA and protein models, and flexible multistate models
 * - user-specified parameter values (e.g. base freqs or alpha parameter)
 * - mixture models (Gamma, FreeRate, P-inv)
 * - Ascertainment bias correction
 * - custom character-to-state mapping
 *
 * For more documentation see: https://github.com/amkozlov/raxml-ng/wiki/Input-data#single-model
 */
class EvolModel
{
public:
  typedef std::unordered_map<int,ParamMode> ParamModeMap;

  EvolModel (DataType data_type = DataType::autodetect, const std::string &model_string = "GTR");
  EvolModel (const std::string &model_string) : EvolModel(DataType::autodetect, model_string) {};

  EvolModel(const EvolModel&) = default;

  /* getters */
  DataType data_type() const { return _data_type; };
  std::string data_type_name() const;
  unsigned int num_states() const { return _num_states; };
  std::string name() const { return _name; };

  /**
   * @brief Mapping between alignment characters and states
   */
  const corax_state_t* charmap() const;

  /**
   * @brief List of *non-ambiguous* state names (eg DNA = A C G T)
   */
  const std::vector<std::string>& state_names() const;

  /**
   * @brief Full map between state id and name, including ambiguities (eg DNA = A C G T M R W S Y K -)
   */
  const StateNameMap& full_state_namemap() const;

  /**
   * @brief Get substitution matrix i
   *
   * For instance, LG4X model has 4 matrices which correspond to 4 evolutionary rates
   * For most other models (LG+G, GTR etc.), there is a single matrix with index 0
   */
  const SubstitutionMatrix submodel(size_t i) const { return _submodels.at(i); };

  /**
   * @brief Rate heterogeneity modes. Possible values:
   * - CORAX_UTIL_MIXTYPE_FIXED  = none
   * - CORAX_UTIL_MIXTYPE_GAMMA  = Gamma model
   * - CORAX_UTIL_MIXTYPE_FREE   = FreeRate model
   */
  unsigned int ratehet_mode() const { return _rate_het; };

  /**
   * @brief Number of rate categories
   */
  unsigned int num_ratecats() const { return _num_ratecats; };

  /**
   * @brief Number of substitution matrices, see submodels() for details
   */
  unsigned int num_submodels() const { return _num_submodels; };

  /**
   * @brief Evolutionary rates for mixture models
   */
  const std::vector<double>& ratecat_rates() const { return _ratecat_rates; };

  /**
   * @brief Mixture weights
   */
  const std::vector<double>& ratecat_weights() const { return _ratecat_weights; };

  /**
   * @brief Mapping between evolutionary rate categories and substitution matrices, e.g.
   *
   * LG+G4 = 0 0 0 0 (all categories use the same matrix)
   * LG4X  = 0 1 2 3 (one individual matrix per category)
   */
  const std::vector<unsigned int>& ratecat_submodels() const { return _ratecat_submodels; };

  /**
   * @brief Gamma estimation method. Possible values:
   * - CORAX_GAMMA_RATES_MEAN
   * - CORAX_GAMMA_RATES_MEDIAN
   */
  int gamma_mode() const { return _gamma_mode; };

  /**
   * @brief Alpha parameter of Gamma model of rate heterogeneity (+G)
   */
  double alpha() const { return _alpha; };

  /**
   * @brief Proportion of invariant sites (+I)
   */
  double pinv() const { return _pinv; };

  /**
   * @brief Per-partition branch length scaler for proportional linkage mode (+B)
   */
  double brlen_scaler() const { return _brlen_scaler; };

  /**
   * @brief Stationary base frequency for i
   */
  const std::vector<double>& base_freqs(unsigned int i) const { return _submodels.at(i).base_freqs(); };
  const std::vector<double>& subst_rates(unsigned int i) const { return _submodels.at(i).subst_rates(); };

  /**
   * @brief Print model definition (and parameters) in compact form, e.g. HKY{1/2.424}+G{1.000}
   *
   * @param print_params whether to include parameter values
   * @param precision number of digits in parameter values
   */
  std::string to_string(bool print_params = false, unsigned int precision = 0) const;

  /**
   * @brief Parameter mask for ML optimization, see CORAX_OPT_PARAM_*
   */
  int params_to_optimize() const;

  /**
   * @brief Get estimation modes for all model parameters
   */
  const ParamModeMap& param_mode() const { return _param_mode; }

  /**
   * @brief Get estimation mode for param, see ParamMode
   */
  ParamMode param_mode(int param) const { return _param_mode.at(param); };

  /**
   * @brief Return true if param is estimated, currently ParamMode::ML || ParamMode::empirical
   */
  bool param_estimated(int param) const;

  /**
   * @brief Get ascertainment bias correction type
   */
  AscBiasCorrection ascbias_type() const { return _ascbias_type; }

  /**
   * @brief Get invariant site counts for Felsenstein or Stamatakis asc. bias correction
   */
  const std::vector<corax_weight_t>& ascbias_weights() const { return _ascbias_weights; }

  /**
   * @brief Size of the CLV entry per alignment site, given in doubles (NOT in bytes)
   *
   * E.g. for DNA data and GTR+G4 model, clv_entry_size = 4*4 = 16 doubles
   */
  size_t clv_entry_size() const { return _num_states * _num_ratecats; }

  /**
   * @brief Number of free parameters in the model, e.g. 5 for GTR or 2 for JC+G
   */
  unsigned int  num_free_params() const;

  /* setters */
  void alpha(double value) { _alpha = value; };
  void pinv(double value) { _pinv = value; };
  void brlen_scaler(double value) { _brlen_scaler = value; };
  void base_freqs(size_t i, const std::vector<double>& value) { _submodels.at(i).base_freqs(value); };
  void subst_rates(size_t i, const std::vector<double>& value) { _submodels.at(i).subst_rates(value); };
  void base_freqs(const std::vector<double>& value) { for (SubstitutionMatrix& s: _submodels) s.base_freqs(value); };
  void subst_rates(const std::vector<double>& value) { for (SubstitutionMatrix& s: _submodels) s.subst_rates(value); };
  void ratecat_rates(std::vector<double> const& value) { _ratecat_rates = value; };
  void ratecat_weights(std::vector<double> const& value) { _ratecat_weights = value; };

  /**
   * @brief Set parameter estimation mode
   */
  void param_mode(int param, ParamMode mode) { _param_mode[param] = mode; };

  /**
   * @brief Set parameter estimation mode only if it is still undefined
   */
  void set_param_mode_default(int param, ParamMode mode)
  {
    if (param_mode(param) == ParamMode::undefined)
      _param_mode[param] = mode;
  };

  /**
   * @brief Parse model from string, e.g. "GTR+G{1.0}"
   */
  void init_from_string(const std::string& model_string);

private:
  std::string _name;
  DataType _data_type;
  unsigned int _num_states;

  std::string _custom_states;
  std::string _custom_gaps;
  bool _custom_case_sensitive;
  std::shared_ptr<corax_state_t> _custom_charmap;
  mutable std::vector<std::string> _state_names;
  mutable StateNameMap _full_state_namemap;

  unsigned int _rate_het;
  unsigned int _num_ratecats;
  unsigned int _num_submodels;
  std::vector<double> _ratecat_rates;
  std::vector<double> _ratecat_weights;
  std::vector<unsigned int> _ratecat_submodels;
  int _gamma_mode;

  double _alpha;
  double _pinv;
  double _brlen_scaler;

  AscBiasCorrection _ascbias_type;
  std::vector<corax_weight_t> _ascbias_weights;

  std::vector<SubstitutionMatrix> _submodels;

  ParamModeMap _param_mode;

  void autodetect_data_type(const std::string& model_name);
  corax_mixture_model_t * init_mix_model(const std::string& model_name);
  void init_model_opts(const std::string& model_opts, const corax_mixture_model_t& mix_model);
  void init_state_names() const;
  void set_user_srates(std::vector<double>& srates, bool normalize = true);
  void set_user_freqs(std::vector<double>& freqs);
};

typedef std::unordered_map<size_t, EvolModel> ModelMap;
typedef std::unordered_map<size_t, EvolModel&> ModelRefMap;
typedef std::unordered_map<size_t, const EvolModel&> ModelCRefMap;

/**
 * @brief Copy model parameter values from corax_partition_t to EvolModel object
 *
 * NOTE: Model definition must be compatible, e.g. number of states, matrices etc.
 */
void assign(EvolModel& model, const corax_partition_t * partition);

/**
 * @brief Copy model parameter values from EvolModel object to corax_partition_t
 *
 * NOTE: Model definition must be compatible, e.g. number of states, matrices etc.
 */
void assign(corax_partition_t * partition, const EvolModel& model);

} // namespace corax
} // namespace model

/**
 * @brief Print "long" model specification (as in raxml-ng log file)
 */
std::ostream& operator<<(std::ostream& stream, const corax::model::EvolModel& m);

#endif /* CORAX_MODEL_EVOLMODEL_HPP_ */
