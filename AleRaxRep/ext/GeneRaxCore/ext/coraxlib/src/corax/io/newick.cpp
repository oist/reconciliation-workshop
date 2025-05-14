#include <cstdlib>
#include <exception>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include "corax/corax.h"
#include "corax/io/newick.hpp"

inline void synchronize_attrs(corax_unode_t *start)
{
  corax_unode_t *cur = start;
  if (start->label)
  {
    while (cur->next != start)
    {
      cur        = cur->next;
      cur->label = start->label;
    }
  }
}

inline void free_label_from_nodes(corax_unode_t *start)
{
  if (start->label)
  {
    free(start->label);
    start->label = nullptr;
  }
  corax_unode_t *cur = start;
  while (cur->next && cur->next != start)
  {
    cur        = cur->next;
    cur->label = nullptr;
  }
}

void delete_unode(corax_unode_t *node)
{
  if (node->label) { free_label_from_nodes(node); }
  free(node);
}

inline size_t close_node_loop(corax_unode_t *start)
{
  size_t         node_count = 1;
  corax_unode_t *cur        = start;
  while (cur->next != nullptr)
  {
    cur = cur->next;
    node_count++;
  }
  cur->next = start;
  return node_count;
}

inline void set_mutual_back_pointers(corax_unode_t *a, corax_unode_t *b)
{
  a->back = b;
  b->back = a;
  if (a->length == 0.0) { a->length = b->length; }
  else { b->length = a->length; }
}

inline void free_node_exception(corax_unode_t *n)
{
  if (n == nullptr) { return; }

  auto s = n->next;
  while (s != n && s != nullptr)
  {
    if (s->back != nullptr) { free_node_exception(s->back); }
    auto tmp = s->next;
    if (s) { delete_unode(s); }
    s = tmp;
  }

  delete_unode(n);
}

inline bool unode_is_rooted(const corax_unode_t *root)
{
  return (root->next && root->next->next == root) ? 1 : 0;
}

static void fill_nodes_recursive(corax_unode_t  *node,
                                 corax_unode_t **array,
                                 unsigned int    array_size,
                                 unsigned int   *tip_index,
                                 unsigned int   *inner_index,
                                 unsigned int    level)
{
  unsigned int index;
  if (!node->next)
  {
    /* tip node */
    index = *tip_index;
    *tip_index += 1;
  }
  else
  {
    /* inner node */
    corax_unode_t *snode = level ? node->next : node;
    do {
      fill_nodes_recursive(
          snode->back, array, array_size, tip_index, inner_index, level + 1);
      snode = snode->next;
    } while (snode != node);

    index = *inner_index;
    *inner_index += 1;
  }

  assert(index < array_size);
  array[index] = node;
}

corax_unode_t *trim_node(corax_unode_t *node)
{
  corax_unode_t *next, *prev;
  next = node->next;
  prev = node;

  while (prev->next != node) { prev = prev->next; }

  prev->next  = next;
  prev->label = node->label;
  node->label = nullptr;
  delete_unode(node);

  return prev->next;
}

corax_unode_t *unode_unroot(corax_unode_t *vroot)
{
  if (vroot->next->next != vroot) { return vroot; }

  corax_unode_t *lchild, *rchild;
  lchild = vroot->back;
  rchild = vroot->next->back;

  double total_len = vroot->length + vroot->next->length;

  lchild->length = total_len;
  rchild->length = total_len;

  lchild->back = rchild;
  rchild->back = lchild;

  delete_unode(vroot->next);
  delete_unode(vroot);

  if (lchild->next != nullptr) { return lchild; }
  return rchild;
}

std::string corax_newick_lexer_t::consume_value_as_string()
{
  std::string tmp;
  std::swap(tmp, _value);
  return tmp;
}

/* WARNING, ALLOCATES MEMORY */
char *corax_newick_lexer_t::consume_value_as_cstring()
{
  char *label = (char *)calloc(
      sizeof(char), (_value.size() + 1) /* Need to include space for the null */
  );

  for (size_t i = 0; i < _value.size(); ++i) { label[i] = _value[i]; }
  _value.clear();
  return label;
}

double corax_newick_lexer_t::consume_value_as_float()
{
  auto   f_str = consume_value_as_string();
  size_t pos   = 0;
  double val   = std::stod(f_str, &pos);
  if (pos != f_str.size())
  {
    throw std::runtime_error{std::string("Float conversion failed around")
                             + describe_position()};
  }
  return val;
}

std::string corax_newick_lexer_t::describe_position() const
{
  std::stringstream builder;
  builder << "position " << _current_index;
  return builder.str();
}

bool corax_newick_lexer_t::is_punct(char c)
{
  return c == '[' || c == ']' || c == '(' || c == ')' || c == ':' || c == ';'
         || c == ',' || c == 0 || c == EOF;
}

std::pair<corax_lexeme_t, size_t> corax_newick_lexer_t::consume_token_pos()
{
  auto start_index = _current_index;
  auto token       = consume();
  return {token, start_index};
}

void corax_newick_lexer_t::skip_whitespace()
{
  while (_current_index < _input.size())
  {
    char c = _input[_current_index];
    if (!std::isspace(c)) { break; }
    _current_index++;
  }
}

void corax_newick_lexer_t::expect(corax_lexeme_t token_type)
{
  auto ret = consume_token_pos();
  if (ret.first != token_type)
  {
    throw std::runtime_error{
        std::string("Got the wrong token type at position ")
        + std::to_string(ret.second + 1) + " was expecting "
        + describe_token(token_type)};
  }
}

std::string corax_newick_lexer_t::describe_token(corax_lexeme_t token_type)
{
  switch (token_type)
  {
  case OPENING_SQUARE_BRACKET:
    return {"opening square bracket"};
  case CLOSING_SQUARE_BRACKET:
    return {"closing square bracket"};
  case OPENING_PAREN:
    return {"opening parenthesis"};
  case CLOSING_PAREN:
    return {"closing parenthesis"};
  case COLON:
    return {"colon"};
  case SEMICOLON:
    return {"semicolon"};
  case COMMA:
    return {"comma"};
  case END:
    return {"end of input"};
  case VALUE:
    return {"either a identifier or a number"};
  default:
    return {"unknown token"};
  }
}

corax_lexeme_t corax_newick_lexer_t::peak()
{
  size_t tmp_index    = _current_index;
  char   current_char = _input[tmp_index++];
  if (is_punct(current_char))
  {
    switch (current_char)
    {
    case '[':
      return OPENING_SQUARE_BRACKET;
    case ']':
      return CLOSING_SQUARE_BRACKET;
    case '(':
      return OPENING_PAREN;
    case ')':
      return CLOSING_PAREN;
    case ':':
      return COLON;
    case ';':
      return SEMICOLON;
    case ',':
      return COMMA;
    case 0:
    case EOF:
      return END;
    default:
      throw std::runtime_error{"The punctuation was unrecognized"};
    }
  }
  else { return VALUE; }
}

corax_lexeme_t corax_newick_lexer_t::consume()
{
  auto token = peak();
  if (token == VALUE)
  {
    // we have a value, so we need to scan until we have found punctuation, or
    // the end of the string
    auto start_itr = _input.begin() + _current_index;
    auto end_itr = start_itr;
    while (char tmp = _input[_current_index])
    {
      if (is_punct(tmp)) { break; }
      end_itr++;
      _current_index++;
    }

    _value = std::string{start_itr, end_itr};
    while (std::isspace(*(_value.end() - 1)))
    {
      _value.resize(_value.size() - 1);
    }
    return token;
  }
  else
  {
    _current_index++;
    skip_whitespace();
    return token;
  }
}

corax_utree_t *corax_newick_parser_t::parse_utree(bool auto_unroot,
                                                  bool allow_rooted)
{
  corax_unode_t *root_node = nullptr;

  try
  {
    root_node = parse_subtree();
    root_node = trim_node(root_node);
    /* We overcounted, because we assume there is an "upper" branch still */
    _edge_count--;
    if (unode_is_rooted(root_node))
    {
      if (!allow_rooted && !auto_unroot)
      {
        throw std::invalid_argument{"Encountered a rooted tree when parsing"};
      }
      if (auto_unroot)
      {
        auto new_root_node = unode_unroot(root_node);
        if (new_root_node != root_node)
        {
          _inner_count--;
          _edge_count--;
          root_node = new_root_node;
        }
      }
    }

    _lexer.expect(SEMICOLON);
    if (!_lexer.at_end())
    {
      throw std::runtime_error{
          "There were extra charcters when we finished parsing"};
    }
  }
  catch (...)
  {
    if (root_node != nullptr)
    {
      free_node_exception(root_node->back);
      free_node_exception(root_node);
    }
    throw;
  }

  auto current_tree         = (corax_utree_t *)calloc(sizeof(corax_utree_t), 1);
  current_tree->tip_count   = _tip_count;
  current_tree->inner_count = _inner_count;
  current_tree->edge_count  = _edge_count;
  current_tree->binary =
      (_inner_count == (_tip_count - (unode_is_rooted(root_node) ? 1 : 2)))
          ? true
          : false;
  size_t node_array_size = _tip_count + _inner_count;
  current_tree->nodes =
      (corax_unode_t **)malloc(sizeof(corax_unode_t *) * node_array_size);

  unsigned int tip_index   = 0;
  unsigned int inner_index = _tip_count;
  fill_nodes_recursive(root_node,
                       current_tree->nodes,
                       node_array_size,
                       &tip_index,
                       &inner_index,
                       /*level=*/0);

  if (tip_index != _tip_count)
  {
    throw std::runtime_error{
        "There was a problem with filling the nodes array"};
  }

  if (inner_index != _tip_count + _inner_count)
  {
    throw std::runtime_error{
        "There was a problem with filling the nodes array"};
  }

  corax_utree_reset_template_indices(root_node, _tip_count);

  current_tree->vroot = root_node;
  return current_tree;
}

corax_unode_t *corax_newick_parser_t::parse_subtree()
{
  auto token = _lexer.peak();
  _edge_count++;

  if (token == OPENING_PAREN)
  {
    auto tmp = parse_internal();
    /* error checking on tmp here */
    return tmp;
  }
  else
  {
    auto tmp = parse_leaf();
    /* error checking on tmp here */
    return tmp;
  }
}

corax_unode_t *corax_newick_parser_t::parse_internal()
{
  corax_unode_t *extra_node   = nullptr;
  corax_unode_t *current_node = nullptr;

  try
  {
    _lexer.expect(OPENING_PAREN);
    extra_node       = (corax_unode_t *)calloc(sizeof(corax_unode_t), 1);
    current_node     = parse_node_set();
    extra_node->next = current_node;
    auto node_count  = close_node_loop(extra_node);
    if (node_count < 3)
    {
      throw std::runtime_error{std::string("Got a singleton subtree around ")
                               + std::string(_lexer.describe_position())};
    }

    _lexer.expect(CLOSING_PAREN);

    parse_node_attrs(extra_node);
    synchronize_attrs(extra_node);

    _inner_count++;
    return extra_node;
  }
  catch (...)
  {
    free_node_exception(extra_node);
    throw;
  }
}

corax_unode_t *corax_newick_parser_t::parse_node_set()
{
  corax_unode_t *current_node = nullptr;
  corax_unode_t *child        = nullptr;
  try
  {
    current_node = (corax_unode_t *)calloc(sizeof(corax_unode_t), 1);
    auto child   = parse_subtree();
    set_mutual_back_pointers(current_node, child);
    auto token = _lexer.peak();
    if (token == COMMA)
    {
      _lexer.consume();
      current_node->next = parse_node_set();
    }
    return current_node;
  }
  catch (...)
  {
    free_node_exception(current_node->back);
    free_node_exception(current_node);
    free_node_exception(child);
    throw;
  }
}

void corax_newick_parser_t::parse_node_attrs(corax_unode_t *current_node)
{
  parse_name(current_node);
  parse_comment();
  parse_length(current_node);
  parse_comment();
}

corax_unode_t *corax_newick_parser_t::parse_leaf()
{
  corax_unode_t *current_node = nullptr;
  try
  {
    current_node = (corax_unode_t *)calloc(sizeof(corax_unode_t), 1);
    parse_node_attrs(current_node);
    if (current_node->label == nullptr)
    {
      throw std::runtime_error{
          std::string("Got a leaf with an empty name around ")
          + std::string(_lexer.describe_position())};
    }
    _tip_count++;
    return current_node;
  }
  catch (...)
  {
    free_node_exception(current_node);
    throw;
  }
}

void corax_newick_parser_t::parse_length(corax_unode_t *current_node)
{
  auto token = _lexer.peak();
  if (token == COLON)
  {
    _lexer.consume();
    _lexer.expect(VALUE);
    current_node->length = parse_number();
  }
}

void corax_newick_parser_t::parse_name(corax_unode_t *current_node)
{
  auto token = _lexer.peak();
  if (token == VALUE)
  {
    _lexer.consume();
    current_node->label = parse_cstring();
  }
}

std::string corax_newick_parser_t::parse_string()
{
  return _lexer.consume_value_as_string();
}

char *corax_newick_parser_t::parse_cstring()
{
  return _lexer.consume_value_as_cstring();
}

double corax_newick_parser_t::parse_number()
{
  return _lexer.consume_value_as_float();
}

void corax_newick_parser_t::parse_comment()
{
  auto token = _lexer.peak();
  if (token == OPENING_SQUARE_BRACKET)
  {
    _lexer.consume();
    _lexer.consume_until(CLOSING_SQUARE_BRACKET);
  }
};

/* Legacy C wrappers */

corax_utree_t *utree_parse_newick_string(std::string newick_string,
                                         bool        auto_unroot,
                                         bool        allow_rooted)
{
  try
  {
    corax_newick_parser_t np(newick_string);
    return np.parse(auto_unroot, allow_rooted);
  }
  catch (std::invalid_argument &e)
  {
    corax_set_error(CORAX_ERROR_INVALID_TREE, e.what());
    return nullptr;
  }
  catch (std::runtime_error &e)
  {
    corax_set_error(CORAX_ERROR_NEWICK_SYNTAX, e.what());
    return nullptr;
  }
}

corax_utree_t *
utree_parse_newick(const char *filename, bool auto_unroot, bool allow_rooted)
{
  std::ifstream newick_file(filename);
  std::string   newick_string((std::istreambuf_iterator<char>(newick_file)),
                            (std::istreambuf_iterator<char>()));
  return utree_parse_newick_string(newick_string, auto_unroot, allow_rooted);
}

CORAX_EXPORT corax_utree_t *corax_utree_parse_newick(const char *filename)
{
  return utree_parse_newick(
      filename, /*auto_unroot=*/false, /*allow_rooted=*/false);
}

CORAX_EXPORT corax_utree_t *
corax_utree_parse_newick_rooted(const char *filename)
{
  return utree_parse_newick(
      filename, /*auto_unroot=*/false, /*allow_rooted=*/true);
}

CORAX_EXPORT corax_utree_t *
corax_utree_parse_newick_unroot(const char *filename)
{
  return utree_parse_newick(
      filename, /*auto_unroot=*/true, /*allow_rooted=*/false);
}

CORAX_EXPORT corax_utree_t *corax_utree_parse_newick_string(const char *s)
{
  return utree_parse_newick_string(std::string(s),
                                   /*auto_unroot=*/false,
                                   /*allow_rooted=*/false);
}

CORAX_EXPORT corax_utree_t *
corax_utree_parse_newick_string_rooted(const char *s)
{
  return utree_parse_newick_string(std::string(s),
                                   /*auto_unroot=*/false,
                                   /*allow_rooted=*/true);
}

CORAX_EXPORT corax_utree_t *
corax_utree_parse_newick_string_unroot(const char *s)
{
  return utree_parse_newick_string(std::string(s),
                                   /*auto_unroot=*/true,
                                   /*allow_rooted=*/false);
}
