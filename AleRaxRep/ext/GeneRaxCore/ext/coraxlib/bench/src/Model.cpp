/*
 * This file was copied from RAxML-NG
 */
#include "Model.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

const std::vector<int> ALL_MODEL_PARAMS = {CORAX_OPT_PARAM_FREQUENCIES, CORAX_OPT_PARAM_SUBST_RATES,
                                      CORAX_OPT_PARAM_PINV, CORAX_OPT_PARAM_ALPHA,
                                      CORAX_OPT_PARAM_FREE_RATES, CORAX_OPT_PARAM_RATE_WEIGHTS,
                                      CORAX_OPT_PARAM_BRANCH_LEN_SCALER,
                                      CORAX_OPT_PARAM_BRANCHES_ITERATIVE};

const std::unordered_map<DataType,unsigned int,EnumClassHash>  DATATYPE_STATES { {DataType::dna, 4},
                                                                            {DataType::protein, 20},
                                                                            {DataType::binary, 2}
                                                                            /*
                                                                            ,

     {DataType::diploid10, 10}
*/

                                                                          };

const std::unordered_map<DataType,const corax_state_t*,EnumClassHash>  DATATYPE_MAPS {
  {DataType::dna, corax_map_nt},
  {DataType::protein, corax_map_aa},
  {DataType::binary, corax_map_bin}
  /*
  ,
  {DataType::diploid10, corax_map_diploid10}
  */
};


// TODO move it out of here
class parse_error {};
static std::string read_option(std::istringstream& s)
{
  std::ostringstream os;
  while (s.peek() != '{' && s.peek() != '+' && s.peek() != EOF)
  {
    os.put(static_cast<char>(toupper(s.get())));
  }
  return os.str();
}

static bool read_param(std::istringstream& s, std::string& val)
{
  if (s.peek() == '{' || s.peek() == '[')
  {
    auto delim = (s.peek() == '{') ? '}' : ']';

    // consume the opening bracket
    s.get();

    char str[50000];
    if (!s.getline(str, 50000, delim))
      throw parse_error();
    val = str;

    return true;
  }
  else
    return false;
}

template<typename T>
static bool read_param(std::istringstream& s, T& val)
{
  if (s.peek() == '{' || s.peek() == '[')
  {
    s.get();
    s >> val;
    auto c = s.get();
    if (c != '}' && c != ']')
      throw parse_error();

    return true;
  }
  else
    return false;
}

template<typename T>
static bool read_param(std::istringstream& s, std::vector<T>& vec)
{
  if (s.peek() == '{' || s.peek() == '[')
  {
    int c = s.get();
    while (s && c != '}' && c != ']')
    {
      T val;
      s >> val;
      vec.push_back(val);
      c = s.get();
    }
    if (c != '}' && c != ']')
      throw parse_error();

    return true;
  }
  else
    return false;
}

template<typename T>
static bool read_param_file(std::istringstream& s, std::vector<T>& vec)
{
  if (s.peek() == '{' || s.peek() == '[')
  {
    // first, try to interpret parameter value as a file name
    auto delim = (s.peek() == '{') ? '}' : ']';
    auto start = s.tellg();

    // consume the opening bracket
    s.get();

    char fname[50000];
    s.getline(fname, 50000, delim);
    std::ifstream fs(fname);
    if (fs.good())
    {
      while (!fs.eof())
      {
        T val;
        fs >> val;
        if (!fs.fail())
          vec.push_back(val);
      }
      return true;
    }
    else
    {
      // if it failed, rewind the stream and try to parse as a list of values
      s.seekg(start);
      return read_param(s, vec);
    }
  }
  return false;
}

template<typename T>
static void print_param(std::ostringstream& s, T val)
{
  s << "{" << val << "}";
}

template<typename T>
static void print_param(std::ostringstream& s, const std::vector<T>& vec)
{
  s << "{";
  size_t i = 0;
  for (auto v: vec)
  {
    s << ((i > 0) ? "/" : "") << v;
    ++i;
  }
  s << "}";
}

Model::Model (DataType data_type, const std::string &model_string) :
    _data_type(data_type), _custom_charmap(nullptr, free)
{
  // RAxML compatibility hack, TODO: generic model name aliases
  const std::string model_string_tmp = model_string == "DNA" ? "GTR+G+F" : model_string;

  init_from_string(model_string_tmp);
}

const corax_state_t * Model::charmap() const
{
  return _custom_charmap ? _custom_charmap.get() : DATATYPE_MAPS.at(_data_type);
}

void Model::init_from_string(const std::string &model_string)
{
  size_t pos = model_string.find_first_of("+{[");
  const std::string model_name = pos == std::string::npos ? model_string : model_string.substr(0, pos);
  const std::string model_opts = pos == std::string::npos ? "" : model_string.substr(pos);

  /* guess data type */
  autodetect_data_type(model_name);

  assert(_data_type != DataType::autodetect);

  /* set number of states based on datatype */
  if (_data_type == DataType::multistate)
  {
    _num_states = corax_util_model_numstates_mult(model_name.c_str());
    _custom_charmap = std::unique_ptr<corax_state_t, void(*)(void*)>(corax_util_model_charmap_mult(_num_states), free);

    //libcorax_check_error("ERROR in model specification |" + model_name + "|");
    assert(_custom_charmap);
  }
  else
    _num_states = DATATYPE_STATES.at(_data_type);

  corax_mixture_model_t * mix_model = init_mix_model(model_name);

  init_model_opts(model_opts, *mix_model);

  corax_util_model_mixture_destroy(mix_model);
}

std::string Model::data_type_name() const
{
  switch (_data_type)
  {
    case DataType::binary:
      return "BIN";
    case DataType::dna:
      return "DNA";
    case DataType::protein:
      return "AA";
    case DataType::diploid10:
      return "GT";
    case DataType::multistate:
      return "MULTI" + std::to_string(_num_states);
    case DataType::autodetect:
      return "AUTO";
    default:
      return "UNKNOWN";
  }
}

void Model::autodetect_data_type(const std::string &model_name)
{
  if (_data_type == DataType::autodetect)
  {
    if (corax_util_model_exists_genotype(model_name.c_str()))
    {
      _data_type = DataType::diploid10;
    }
    else if (corax_util_model_exists_mult(model_name.c_str()))
    {
      _data_type = DataType::multistate;
    }
    else if (model_name == "BIN")
    {
      _data_type = DataType::binary;
    }
    else if (corax_util_model_exists_protein(model_name.c_str()) ||
             corax_util_model_exists_protmix(model_name.c_str()))
    {
      _data_type = DataType::protein;
    }
    else
    {
      _data_type = DataType::dna;
    }
  }
}

corax_mixture_model_t * Model::init_mix_model(const std::string &model_name)
{
  const char * model_cstr = model_name.c_str();
  corax_mixture_model_t * mix_model = nullptr;

  if (corax_util_model_exists_protmix(model_cstr))
  {
    mix_model = corax_util_model_info_protmix(model_cstr);
  }
  else
  {
    corax_subst_model_t * modinfo =  NULL;

    /* initialize parameters from the model */
    if (_data_type == DataType::protein)
    {
      modinfo =  corax_util_model_info_protein(model_cstr);
    }
    else if (_data_type == DataType::dna)
    {
      modinfo =  corax_util_model_info_dna(model_cstr);
    }
    else if (_data_type == DataType::binary)
    {
      modinfo =  corax_util_model_create_custom("BIN", 2, NULL, NULL, NULL, NULL);
    }
    else if (_data_type == DataType::diploid10)
    {
      modinfo =  corax_util_model_info_genotype(model_cstr);
    }
    else if (_data_type == DataType::multistate)
    {
      modinfo =  corax_util_model_info_mult(model_cstr);
    }

    // TODO: user models must be defined explicitly
//    /* pre-defined model not found; assume model std::string encodes rate symmetries */
//    if (!modinfo)
//      modinfo =  corax_util_model_create_custom("USER", _num_states, NULL, NULL, model_cstr, NULL);

    if (!modinfo)
    {
      /*
      if (corax_errno)
        libcorax_check_error("ERROR model initialization |" + model_name + "|");
      else*/
        throw std::runtime_error("Invalid model name: " + model_name);
    }

    /* create pseudo-mixture with 1 component */
    mix_model = corax_util_model_mixture_create(modinfo->name, 1, &modinfo, NULL, NULL,
                                                  CORAX_UTIL_MIXTYPE_FIXED);

    corax_util_model_destroy(modinfo);
  }

  return mix_model;
}

void Model::init_model_opts(const std::string &model_opts, const corax_mixture_model_t& mix_model)
{
  _gamma_mode = CORAX_GAMMA_RATES_MEAN;
  _alpha = 1.0;
  _pinv = 0.0;
  _brlen_scaler = 1.0;
  _name = std::string(mix_model.name);

  _ascbias_type = AscBiasCorrection::none;
  _ascbias_weights.clear();

  _ratecat_rates.clear();
  _ratecat_weights.clear();

  /* set rate heterogeneity defaults from model */
  _num_ratecats = mix_model.ncomp;
  _num_submodels = mix_model.ncomp;
  _rate_het = static_cast<unsigned int>(mix_model.mix_type);

  /* allocate space for all subst matrices */
  for (size_t i = 0; i < mix_model.ncomp; ++i)
    _submodels.emplace_back(*mix_model.models[i]);

  /* set default param optimization modes */
  for (auto param: ALL_MODEL_PARAMS)
    _param_mode[param] = ParamValue::undefined;

  _param_mode[CORAX_OPT_PARAM_FREQUENCIES] =
      mix_model.models[0]->freqs ? ParamValue::model : ParamValue::ML;

  _param_mode[CORAX_OPT_PARAM_SUBST_RATES] =
      mix_model.models[0]->rates ? ParamValue::model : ParamValue::ML;

  const char *s = model_opts.c_str();

  std::istringstream ss(model_opts);

  try
  {
    doubleVector user_srates;
    if (read_param_file(ss, user_srates))
    {
      // TODO support multi-matrix models
      if (_submodels.size() > 0)
        std::runtime_error("User-defined rates for multi-matrix models are not supported yet!");

      auto smodel = _submodels[0];
      auto num_uniq_rates = smodel.num_uniq_rates();
      if (user_srates.size() != num_uniq_rates)
      {
        throw std::runtime_error("Invalid number of substitution rates specified: " +
                            std::to_string(user_srates.size()) + " (expected: " +
                            std::to_string(num_uniq_rates) + ")\n");
      }

      // normalize the rates
      auto last_rate = smodel.rate_sym().empty() ?
        user_srates.back() : user_srates[static_cast<unsigned int>(smodel.rate_sym().back())];
      for (auto& r: user_srates)
        r /= last_rate;

      for (auto& m: _submodels)
        m.uniq_subst_rates(user_srates);

      _param_mode[CORAX_OPT_PARAM_SUBST_RATES] = ParamValue::user;
    }
  }
  catch(parse_error& e)
  {
    throw std::runtime_error(std::string("Invalid substitution rate specification: ") + s);
  }

  // skip "+"
  ss.get();

  /* parse the rest and set additional model params */
  ParamValue param_mode;
  while(ss.peek() != EOF)
  {
    // set current std::string position, for error reporting
    s = model_opts.c_str() + ss.tellg();

    switch(toupper(ss.get()))
    {
      case EOF:
      case '\0':
        // end of model options
        break;
      case '+':
        // proceed to the next token
        continue;
      case 'A':
      {
        try
        {
          const std::string asc_str = "A" + read_option(ss);
          if (asc_str == "ASC_LEWIS")
          {
            _ascbias_type = AscBiasCorrection::lewis;
          }
          else if (asc_str == "ASC_FELS")
          {
            _ascbias_type = AscBiasCorrection::felsenstein;
            WeightType w;
            if (read_param(ss, w))
            {
              _ascbias_weights.resize(_num_states, 0);
              _ascbias_weights[0] = w;
            }
            else
              throw parse_error();
          }
          else if (asc_str == "ASC_STAM")
          {
            _ascbias_type = AscBiasCorrection::stamatakis;
            WeightVector v;
            if (read_param(ss, v))
            {
              if (v.size() == _num_states)
                _ascbias_weights.assign(v.cbegin(), v.cend());
              else
                throw parse_error();
            }
            else
              throw parse_error();
          }
          else
            throw parse_error();
        }
        catch (parse_error& e)
        {
          throw std::runtime_error(std::string("Invalid ascertainment bias correction specification: ") + s);
        }
        break;
      }
      case 'B':
        try
        {
          int mode = ss.peek() == '{' ? 'U' : toupper(ss.get());
          switch (mode)
          {
            case EOF:
            case '\0':
            case '+':
            case 'O':
              param_mode = ParamValue::ML;
              break;
            case 'U':
              if (read_param(ss, _brlen_scaler))
                param_mode = ParamValue::user;
              else
                throw parse_error();
              break;
            default:
              throw parse_error();
          }
          _param_mode[CORAX_OPT_PARAM_BRANCH_LEN_SCALER] = param_mode;
        }
        catch(parse_error& e)
        {
          throw std::runtime_error(std::string("Invalid branch length scaler specification: ") + s);
        }
        break;
      case 'F':
        try
        {
          int mode = ss.peek() == '{' ? 'U' : toupper(ss.get());
          switch (mode)
          {
            case EOF:
            case '\0':
            case '+':
            case 'C':
              param_mode = ParamValue::empirical;
              break;
            case 'O':
              param_mode = ParamValue::ML;
              break;
            case 'E':
              param_mode = ParamValue::equal;
              break;
            case 'U':
            {
              param_mode = ParamValue::user;

              /* for now, read single set of frequencies */
              doubleVector user_freqs;
              if (read_param_file(ss, user_freqs))
              {
                if (user_freqs.size() != _num_states)
                {
                  throw std::runtime_error("Invalid number of user frequencies specified: " +
                           std::to_string(user_freqs.size()) + "\n" +
                           "Number of frequencies must be equal to the number of states: " +
                           std::to_string(_num_states) + "\n");
                }

                double sum = 0.;
                bool invalid = false;
                for (auto v: user_freqs)
                {
                  invalid |= (v <= 0. || v >= 1.);
                  sum += v;
                }

                if (invalid)
                {
                  throw std::runtime_error("Invalid base frequencies specified! "
                      "Frequencies must be positive numbers between 0. and 1.");
                }

                // normalize freqs
                for (auto& f: user_freqs)
                  f /= sum;

                for (auto& m: _submodels)
                  m.base_freqs(user_freqs);
              }
              else
                throw parse_error();
            }
            break;
            default:
              throw parse_error();
          }
          _param_mode[CORAX_OPT_PARAM_FREQUENCIES] = param_mode;
        }
        catch(parse_error& e)
        {
          throw std::runtime_error(std::string("Invalid frequencies specification: ") + s);
        }
        break;
      case 'I':
        try
        {
          int mode = ss.peek() == '{' ? 'U' : toupper(ss.get());
          switch (mode)
          {
            case EOF:
            case '\0':
            case '+':
            case 'O':
              param_mode = ParamValue::ML;
              break;
            case 'C':
              param_mode = ParamValue::empirical;
              break;
            case 'U':
              if (read_param(ss, _pinv))
                param_mode = ParamValue::user;
              else
                throw parse_error();
              break;
            default:
              throw parse_error();
          }
          _param_mode[CORAX_OPT_PARAM_PINV] = param_mode;
        }
        catch(parse_error& e)
        {
          throw std::runtime_error(std::string("Invalid p-inv specification: ") + s);
        }
        break;
      case 'G':
        try
        {
          /* allow to override mixture ratehet mode for now */
          _rate_het = CORAX_UTIL_MIXTYPE_GAMMA;
          if (isdigit(ss.peek()))
          {
            ss >> _num_ratecats;
          }
          else if (_num_ratecats == 1)
            _num_ratecats = 4;

          if (ss.peek() == 'a' || ss.peek() == 'A')
          {
            ss.get();
            _gamma_mode = CORAX_GAMMA_RATES_MEDIAN;
          }
          else if (ss.peek() == 'm' || ss.peek() == 'M')
          {
            ss.get();
            _gamma_mode = CORAX_GAMMA_RATES_MEAN;
          }

          if (read_param(ss, _alpha))
          {
            _param_mode[CORAX_OPT_PARAM_ALPHA] = ParamValue::user;
          }
        }
        catch(parse_error& e)
        {
          throw std::runtime_error(std::string("Invalid GAMMA specification: ") + s);
        }
        break;
      case 'M':
        try
        {
           std::string state_chars, gap_chars;
           int case_sensitive = 1;

           if (tolower(ss.peek()) == 'i')
           {
             ss.get();
             case_sensitive = 0;
           }

           if (!read_param(ss, state_chars))
             throw parse_error();

           read_param(ss, gap_chars);

           _custom_charmap = std::unique_ptr<corax_state_t, void(*)(void*)>(
               corax_util_charmap_create(_num_states,
                                          state_chars.c_str(),
                                          gap_chars.c_str(),
                                          case_sensitive),
               free);
           if (!_custom_charmap)
           {
             assert(corax_errno);
             throw parse_error();
           }
        }
        catch(parse_error& e)
        {
          throw std::runtime_error(std::string("Invalid character std::map specification: ") + s);
        }
        break;
      case 'R':
        _rate_het = CORAX_UTIL_MIXTYPE_FREE;
        if (isdigit(ss.peek()))
        {
          ss >> _num_ratecats;
        }
        else if (_num_ratecats == 1)
          _num_ratecats = 4;

        try
        {
          doubleVector v;
          if (read_param(ss, v))
          {
            if (v.size() != _num_ratecats)
            {
              throw std::runtime_error("Invalid number of free rates specified: " +
                                  std::to_string(v.size()) + " (expected: " +
                                  std::to_string(_num_ratecats) + ")\n");
            }

            // TODO: maybe allow to optimize rates and weights separately
            _param_mode[CORAX_OPT_PARAM_RATE_WEIGHTS] = ParamValue::user;
            _param_mode[CORAX_OPT_PARAM_FREE_RATES] = ParamValue::user;

            _ratecat_rates = v;

            v.clear();
            if (read_param(ss, v))
            {
              if (v.size() != _num_ratecats)
              {
                throw std::runtime_error("Invalid number of rate weights specified: " +
                                    std::to_string(v.size()) + " (expected: " +
                                    std::to_string(_num_ratecats) + ")\n");
              }


              // normalize weights
              double sum = 0;
              for (auto w: v)
                sum += w;

              for (auto& w: v)
                w /= sum;

              _ratecat_weights = v;
            }
            else
              _ratecat_weights.assign(_num_ratecats, 1.0 / _num_ratecats);

            // normalize weights + rates
            double sum_weightrates = 0.0;
            for (size_t i = 0; i < _num_ratecats; ++i)
              sum_weightrates += _ratecat_rates[i] * _ratecat_weights[i];

            for (auto& r: _ratecat_rates)
              r /= sum_weightrates;
          }
        }
        catch(parse_error& e)
        {
          throw std::runtime_error(std::string("Invalid FreeRate specification: ") + s);
        }
        break;
      default:
        throw std::runtime_error("Wrong model specification: " + model_opts);
    }
  }

  switch (_param_mode.at(CORAX_OPT_PARAM_FREQUENCIES))
  {
    case ParamValue::user:
    case ParamValue::empirical:
    case ParamValue::model:
      /* nothing to do here */
      break;
    case ParamValue::equal:
    case ParamValue::ML:
      /* use equal frequencies as s a starting value for ML optimization */
      for (auto& m: _submodels)
        m.base_freqs(doubleVector(_num_states, 1.0 / _num_states));
      break;
    default:
      assert(0);
  }

  switch (_param_mode.at(CORAX_OPT_PARAM_SUBST_RATES))
  {
    case ParamValue::user:
      break;
    case ParamValue::empirical:
    case ParamValue::model:
      /* nothing to do here */
      break;
    case ParamValue::equal:
    case ParamValue::ML:
      /* use equal rates as s a starting value for ML optimization */
      for (auto& m: _submodels)
        m.subst_rates(doubleVector(m.num_rates(), 1.0));
      break;
    default:
      assert(0);
  }

  /* default: equal rates & weights */
  if (_ratecat_rates.empty())
    _ratecat_rates.assign(_num_ratecats, 1.0);
  if (_ratecat_weights.empty())
    _ratecat_weights.assign(_num_ratecats, 1.0 / _num_ratecats);
  _ratecat_submodels.assign(_num_ratecats, 0);

  if (_num_ratecats > 1)
  {
    /* init rate & weights according to the selected mode */
    switch (_rate_het)
    {
      case CORAX_UTIL_MIXTYPE_FIXED:
        assert(_num_ratecats == mix_model.ncomp);
        /* set rates and weights from the mixture model definition */
        _ratecat_rates.assign(mix_model.mix_rates, mix_model.mix_rates + _num_ratecats);
        _ratecat_weights.assign(mix_model.mix_weights, mix_model.mix_weights + _num_ratecats);
        break;

      case CORAX_UTIL_MIXTYPE_GAMMA:
        /* compute the discretized category rates from a gamma distribution
           with given alpha shape and store them in rate_cats  */
        corax_compute_gamma_cats(_alpha, _num_ratecats, _ratecat_rates.data(), _gamma_mode);
        if (_param_mode[CORAX_OPT_PARAM_ALPHA] == ParamValue::undefined)
          _param_mode[CORAX_OPT_PARAM_ALPHA] = ParamValue::ML;
        break;

      case CORAX_UTIL_MIXTYPE_FREE:
        if (_param_mode[CORAX_OPT_PARAM_FREE_RATES] == ParamValue::undefined)
        {
          /* use GAMMA rates as initial values -> can be changed */
          corax_compute_gamma_cats(_alpha, _num_ratecats, _ratecat_rates.data(), _gamma_mode);
          _param_mode[CORAX_OPT_PARAM_FREE_RATES] = ParamValue::ML;
        }
        if (_param_mode[CORAX_OPT_PARAM_RATE_WEIGHTS] == ParamValue::undefined)
          _param_mode[CORAX_OPT_PARAM_RATE_WEIGHTS] = ParamValue::ML;
        break;

      default:
        throw std::runtime_error("Unknown rate heterogeneity model");
    }

    /* link rate categories to corresponding mixture components (R-matrix + freqs)*/
    if (_num_submodels == _num_ratecats)
    {
      for (unsigned int i = 0; i < _num_ratecats; ++i)
        _ratecat_submodels[i] = i;
    }
  }
}

std::string Model::to_string(bool print_params, unsigned int precision) const
{
  std::ostringstream model_string;
  model_string << name();

  auto out_param_mode = _param_mode;
  if (print_params)
  {
    for (auto& entry: out_param_mode)
      entry.second = (entry.second == ParamValue::ML) ? ParamValue::user : entry.second;
  }

  if (precision)
    model_string << std::fixed << std::setprecision(static_cast<int>(precision));

  if (out_param_mode.at(CORAX_OPT_PARAM_SUBST_RATES) == ParamValue::user)
    print_param(model_string, submodel(0).uniq_subst_rates());

  switch(out_param_mode.at(CORAX_OPT_PARAM_FREQUENCIES))
  {
    case ParamValue::empirical:
      model_string << "+FC";
      break;
    case ParamValue::ML:
      model_string << "+FO";
      break;
    case ParamValue::equal:
      model_string << "+FE";
      break;
    case ParamValue::user:
    {
      model_string << "+FU";
      print_param(model_string, base_freqs(0));
    }
    break;
    default:
      break;
  }

  switch(out_param_mode.at(CORAX_OPT_PARAM_PINV))
  {
    case ParamValue::empirical:
      model_string << "+IC";
      break;
    case ParamValue::ML:
      model_string << "+I";
      break;
    case ParamValue::user:
      model_string << "+IU{" << _pinv << "}";
      break;
    default:
      break;
  }

  if (_num_ratecats > 1)
  {
    if (_rate_het == CORAX_UTIL_MIXTYPE_GAMMA)
    {
      model_string << "+G" << _num_ratecats;
      model_string << (_gamma_mode == CORAX_GAMMA_RATES_MEDIAN ? "a" : "m");
      if (out_param_mode.at(CORAX_OPT_PARAM_ALPHA) == ParamValue::user)
        model_string << "{" << _alpha << "}";
    }
    else if (_rate_het == CORAX_UTIL_MIXTYPE_FREE)
    {
      model_string << "+R" << _num_ratecats;
      if (out_param_mode.at(CORAX_OPT_PARAM_FREE_RATES) == ParamValue::user)
      {
        print_param(model_string, _ratecat_rates);
        print_param(model_string, _ratecat_weights);
      }
    }
  }

  switch(out_param_mode.at(CORAX_OPT_PARAM_BRANCH_LEN_SCALER))
  {
    case ParamValue::ML:
      model_string << "+B";
      break;
    case ParamValue::user:
      model_string << "+BU{" << _brlen_scaler << "}";
      break;
    default:
      break;
  }

  switch(_ascbias_type)
  {
    case AscBiasCorrection::lewis:
      model_string << "+ASC_LEWIS";
      break;
    case AscBiasCorrection::felsenstein:
      model_string << "+ASC_FELS";
      print_param(model_string, _ascbias_weights.at(0));
      break;
    case AscBiasCorrection::stamatakis:
      model_string << "+ASC_STAM";
      print_param(model_string, _ascbias_weights);
      break;
    default:
      break;
  }

  return model_string.str();
}

int Model::params_to_optimize() const
{
  int params_to_optimize = 0;

  for (auto param: ALL_MODEL_PARAMS)
  {
    if (_param_mode.at(param) == ParamValue::ML)
      params_to_optimize |= param;
  }

  return params_to_optimize;
}

void assign(Model& model, const corax_partition_t * partition)
{
  if (model.num_states() == partition->states &&
      model.num_submodels() == partition->rate_matrices)
  {
    model.pinv(partition->prop_invar[0]);
    for (size_t i = 0; i < model.num_submodels(); ++i)
    {
      model.base_freqs(i, doubleVector(partition->frequencies[i],
                                       partition->frequencies[i] + partition->states));

      size_t n_subst_rates = CORAX_SUBST_RATE_COUNT(partition->states);
      model.subst_rates(i, doubleVector(partition->subst_params[i],
                                        partition->subst_params[i] + n_subst_rates));
    }

    if (partition->rate_cats > 1)
    {
      model.ratecat_rates(doubleVector(partition->rates,
                                       partition->rates + partition->rate_cats));
      model.ratecat_weights(doubleVector(partition->rate_weights,
                                         partition->rate_weights + partition->rate_cats));
    }
  }
  else
    throw std::runtime_error("incompatible partition!");
}

void assign(corax_partition_t * partition, const Model& model)
{
  if (model.num_states() == partition->states &&
      model.num_submodels() == partition->rate_matrices)
  {
    /* set rate categories & weights */
    corax_set_category_rates(partition, model.ratecat_rates().data());
    corax_set_category_weights(partition, model.ratecat_weights().data());

    /* now iterate over rate matrices and set all params */
    for (unsigned int i = 0; i < partition->rate_matrices; ++i)
    {
      /* set base frequencies */
      assert(!model.base_freqs(i).empty());
      corax_set_frequencies(partition, i, model.base_freqs(i).data());

      /* set substitution rates */
      assert(!model.subst_rates(i).empty());
      corax_set_subst_params(partition, i, model.subst_rates(i).data());

      /* set p-inv value */
      corax_update_invariant_sites_proportion (partition, i, model.pinv());
    }
  }
  else
    throw std::runtime_error("incompatible partition!");
}
