#ifndef CORAX_IO_NEWICK_HPP
#define CORAX_IO_NEWICK_HPP

#include <iostream>
#include "corax/corax.h"

/** @defgroup corax_newick_parser_t
 * Module containing the newick parser and lexer
 */

/**
 * Lexeme types for the newick parser.
 * @ingroup corax_newick_parser_t
 */
enum corax_lexeme_t
{
  OPENING_SQUARE_BRACKET,
  CLOSING_SQUARE_BRACKET,
  OPENING_PAREN,
  CLOSING_PAREN,
  COLON,
  SEMICOLON,
  COMMA,
  VALUE,
  END
};

/**
 * A lexer for the newick parser. In general, users should not consume this
 * directly, but instead use the `corax_newick_parser_t`
 * @ingroup corax_newick_parser_t
 */
class corax_newick_lexer_t
{
public:
  corax_newick_lexer_t(std::string input) :
      _input{std::move(input)}, _current_index{0} {};

  corax_lexeme_t consume();
  corax_lexeme_t peak();

  std::string consume_value_as_string();

  /* WARNING, ALLOCATES MEMORY */
  /**
   * WARNING, ALLOCATES MEMORY.
   *
   * The user of this function is responsible for deallocating the memory
   * allocated by this function
   */
  char *consume_value_as_cstring();

  double consume_value_as_float();

  std::string describe_position() const;

  void expect(corax_lexeme_t token_type);

  void consume_until(corax_lexeme_t token_type)
  {
    while (token_type != consume()) {}
  }

  bool at_end() { return _input.size() == _current_index; }

private:
  bool is_punct(char c);

  std::pair<corax_lexeme_t, size_t> consume_token_pos();

  std::string describe_token(corax_lexeme_t token_type);

  void skip_whitespace();

  std::string _input;
  std::string _value;
  size_t      _current_index;
};

/**
 * Parser for a newick encoded tree.
 *
 * Using the following grammar:
 * ```
 * <tree> ::=
 *     <subtree> ";"
 * <subtree> ::=
 *     <leaf> |
 *     <internal>
 * <internal> ::=
 *     "(" <node_set> ")" <node_attrs>
 * <node_set> ::=
 *     <node> |
 *     <node> "," <node_set>
 * <node> ::=
 *     <subtree> <node_attrs>
 * <node_attrs> ::=
 *     <name> <length> <comment>
 * <leaf> ::=
 *     <node_attrs>
 * <length> ::=
 *     ":" <number> |
 *     <empty>
 * <name> ::=
 *     <string> |
 *     <empty>
 * <string> ::=
 *     anything but punctuation
 * <number> ::=
 *     [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?
 * <comment> ::=
 *     "[" .* "]" |
 *     <empty>
 * <empty> ::=
 *     ""
 * ```
 *
 * @ingroup corax_newick_parser_t
 */
class corax_newick_parser_t
{
public:
  /**
   * @param input Newick string to be parsed
   */
  corax_newick_parser_t(std::string input) :
      _lexer{std::move(input)},
      _tip_count{0},
      _inner_count{0},
      _edge_count{0} {};

  /**
   * Parse the newick tree which was passed to the constructor.
   */
  corax_utree_t *parse() { return parse_utree(false, false); }

  /**
   * Parse the newick tree which was passed to the constructor. Includes flags
   * for unrooting the tree, and throwing an error when the tree is rooted.
   *
   * @param auto_unroot Unroot the tree if a rooted tree is encountered. Implies
   * allow_rooted.
   * @param allow_rooted Don't throw an error when encountering a rooted tree.
   */
  corax_utree_t *parse(bool auto_unroot, bool allow_rooted)
  {
    return parse_utree(auto_unroot, allow_rooted);
  }

private:
  /**
   * Function corresponding to the rule
   *
   * ```
   * <tree> ::=
   *     <subtree> ";"
   * ```
   */
  corax_utree_t *parse_utree(bool auto_unroot, bool allow_rooted);

  /**
   * Function corresponding to the rule
   *
   * ```
   * <subtree> ::=
   *     <leaf> |
   *     <internal>
   * ```
   */
  corax_unode_t *parse_subtree();

  /**
   * Function corresponding to the rule
   *
   * ```
   * <internal> ::=
   *     "(" <node_set> ")" <node_attrs>
   *
   * Allocates memory by creating a `corax_unode_t`
   * ```
   */
  corax_unode_t *parse_internal(); // creates node

  /**
   * Function corresponding to the rules
   * ```
   * <node_set> ::=
   *     <node> |
   *     <node> "," <node_set>
   * <node> ::=
   *     <subtree> <node_attrs>
   * ```
   */
  corax_unode_t *parse_node_set();

  /**
   * Function corresponding to the rule
   * ```
   * <node_attrs> ::=
   *     <name> <length> <comment>
   * ```
   */
  void parse_node_attrs(corax_unode_t *current_node);

  /**
   * Function corresponding to the rule
   * ```
   * <leaf> ::=
   *     <node_attrs>
   * ```
   *
   * Allocates memory by creating a `corax_unode_t`
   */
  corax_unode_t *parse_leaf(); // creates node

  /**
   * Function corresponding to the rule
   * ```
   * <length> ::=
   *     ":" <number> |
   *     <empty>
   * ```
   */
  void parse_length(corax_unode_t *current_node);

  /**
   * Function corresponding to the rule
   * ```
   * <name> ::=
   *     <string> |
   *     <empty>
   * ```
   */
  void parse_name(corax_unode_t *current_node);

  /**
   * Function corresponding to the rule
   * ```
   * <string> ::=
   *     anything but punctuation
   * ```
   */
  std::string parse_string();

  /**
   * Function corresponding to the rule
   * ```
   * <string> ::=
   *     anything but punctuation
   * ```
   */
  char *parse_cstring();

  /**
   * Function corresponding to the rule
   * ```
   * <number> ::=
   *     [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?
   * ```
   */
  double parse_number();

  /**
   * Function corresponding to the rule
   * ```
   * <comment> ::=
   *     "[" .* "]" |
   *     <empty>
   * ```
   */
  void parse_comment();

  /* member variables */
  corax_newick_lexer_t _lexer;
  size_t               _tip_count;
  size_t               _inner_count;
  size_t               _edge_count;
};
#endif // CORAX_IO_NEWICK_HPP
