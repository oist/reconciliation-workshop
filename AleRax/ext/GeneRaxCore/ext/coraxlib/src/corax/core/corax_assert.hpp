#pragma once
#ifndef _CORAX_ASSERT_HPP

#include <sstream>
#include <iostream>
#include <string>

// The exception throws by corax_always_assert() when
// CORAX_ASSERT_USE_EXCEPTIONS is defined.
#ifdef CORAX_ASSERT_USE_EXCEPTIONS
class corax_assert_failed_exception : public std::exception
{
public:
  corax_assert_failed_exception(char const *const message) : _message(message)
  {
  }
  corax_assert_failed_exception(const std::string message) : _message(message)
  {
  }

  const char *what() const noexcept override { return _message.c_str(); }

private:
  const std::string _message;
};
#endif

// A function to build the error message of failed assertions.
inline std::string corax_assert_message(const char *expression,
                                        const char *message,
                                        const char *file,
                                        int         line)
{
  std::stringstream ss;
  ss << "Assertion '" << expression << "' failed in " << file << " line "
     << line;
  if (message != nullptr && message[0] != '\0') { ss << ": " << message; }
  return ss.str();
}

#ifndef NDEBUG
// Expanding the macro into a `do { ... } while(false)` pseudo-loop makes it
// "act like a statement". Without this, e.g. a macro inside an if-statement
// without surrounding it in braces would have unexpected effects.

#define CORAX_ASSERT_IMPL(Expression, Message)                                 \
  do {                                                                         \
    if (!(Expression))                                                         \
    {                                                                          \
      std::cerr << corax_assert_message(                                       \
          #Expression, Message, __FILE__, __LINE__)                            \
                << std::endl;                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (false)
#else
#define CORAX_ASSERT_IMPL(Expression, Message)
#endif

#ifndef CORAX_ASSERT_USE_EXCEPTIONS
#define CORAX_ALWAYS_ASSERT_IMPL(Expression, Message)                          \
  do {                                                                         \
    if (!(Expression))                                                         \
    {                                                                          \
      std::cerr << corax_assert_message(                                       \
          #Expression, Message, __FILE__, __LINE__)                            \
                << std::endl;                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (false)
#else
#define CORAX_ALWAYS_ASSERT_IMPL(Expression, Message)                          \
  do {                                                                         \
    if (!(Expression))                                                         \
    {                                                                          \
      throw corax_assert_failed_exception(                                     \
          corax_assert_message(#Expression, Message, __FILE__, __LINE__));     \
    }                                                                          \
  } while (false)
#endif

// We want to support an additional message for corax_assert and thus need an
// overload with an additional parameter. corax_always_assert() always requires
// an error messsage.
#define CORAX_ASSERT_2_ARGS(Expression, Message)                               \
  CORAX_ASSERT_IMPL(Expression, Message)
#define CORAX_ASSERT_1_ARGS(Expression) CORAX_ASSERT_2_ARGS(Expression, nullptr)
#define CORAX_ASSERT_GET_3rd_ARG(Arg1, Arg2, Arg3, ...) Arg3
#define CORAX_ASSERT_OVERLOAD_CHOOSER(...)                                     \
  CORAX_ASSERT_GET_3rd_ARG(                                                    \
      __VA_ARGS__, CORAX_ASSERT_2_ARGS, CORAX_ASSERT_1_ARGS, DUMMY)
#define corax_assert(...)                                                      \
  CORAX_ASSERT_OVERLOAD_CHOOSER(__VA_ARGS__)(__VA_ARGS__)
#define corax_always_assert(Expression, Message)                               \
  CORAX_ALWAYS_ASSERT_IMPL(Expression, Message)

#endif // _CORAX_ASSERT_HPP
