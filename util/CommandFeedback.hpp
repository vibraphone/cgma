#ifndef CUBIT_COMMANDFEEDBACK_HPP
#define CUBIT_COMMANDFEEDBACK_HPP

#include <typeinfo>
#include <cassert>

// Base class for all command feedback.
class CommandFeedback
{
public:
  typedef const std::type_info& Type;

  CommandFeedback()
  {}

  virtual ~CommandFeedback()
  {}

  virtual Type this_type() const = 0;

  template <typename target_type>
    target_type& get_reference()
  {
    assert(this->this_type() == target_type::type());
    return *(dynamic_cast<target_type*>(this));
  }

  template <typename target_type>
    const target_type& get_reference() const
  {
    assert(this->this_type() == target_type::type());
    return *(dynamic_cast<const target_type*>(this));
  }
};

// The rest of this file is an example to show how to use CommandFeedback
#include "CubitString.hpp"

class DeprecatedCommandError : public CommandFeedback
{
public:
  DeprecatedCommandError(const CubitString &some_string)
    : myInt(27), myString(some_string) 
  {}

  int wow_an_integer() const
  { return myInt; }
  CubitString wow_a_string() const
  { return myString; }

  static CommandFeedback::Type type()
  { return typeid(DeprecatedCommandError); }

  CommandFeedback::Type this_type() const
  { return DeprecatedCommandError::type(); }

private:
  int myInt;
  CubitString myString;
};

#endif
