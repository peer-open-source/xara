//
// Description: Definition of the Hash class, a compile-time string hash
// function. This is useful for automatically generating unique
// class tags at compile time.
//
// Inspired by the stackoverflow answers here: 
// - https://stackoverflow.com/a/23684632
// - https://stackoverflow.com/q/48896142
//
// EXAMPLE:
//
//     using namespace OpenSees::Hash::literals;
//     void one() {} void two() {} void other() {}
//
//     void foo( const std::string& value )
//     {
//       switch( Hash::hash(value) )
//       {
//         case "one"_hash: one(); break;
//         case "two"_hash: two(); break;
//         // many more cases
//         default: other(); break;
//       }
//     }
//
//
// Claudio Perez
// Spring 2023
//
#ifndef OpenSeesHash_H
#define OpenSeesHash_H
#include <string>
#include <utility>
#include <type_traits>


namespace OpenSees {
namespace Hash {
  typedef std::size_t hash_t;
  // typedef hash_t hash_t;
  //
  template<class>struct hasher;
  template<>
  struct hasher<std::string> {
    hash_t constexpr operator()(char const *input)const {
      return *input ?
        static_cast<unsigned int>(*input) + 33 * (*this)(input + 1) :
        5381;
    }
    hash_t operator()( const std::string& str ) const {
      return (*this)(str.c_str());
    }
  };
  template<typename T>
  hash_t constexpr hash(T&& t) {
    return hasher< typename std::decay<T>::type >()(std::forward<T>(t));
  }
  inline namespace literals {
    hash_t constexpr operator "" _hash(const char* s, size_t) {
      return hasher<std::string>()(s);
    }
  }
}
}
#endif // OpenSeesHash_H
