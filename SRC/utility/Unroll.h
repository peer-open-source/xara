//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//----------------------------------------------------------------------------//
#pragma once
#include <utility>
#include <cstddef>

#if defined(__GNUC__) || defined(__clang__)
#  define XARA_UNROLL_INLINE [[gnu::always_inline]] inline
#elif defined(_MSC_VER)
#  define XARA_UNROLL_INLINE __forceinline
#else
#  define XARA_UNROLL_INLINE inline
#endif

namespace OpenSees {
namespace
{
  template<std::size_t N>
  struct num { static constexpr std::size_t value = N; };

  // Helper that unpacks the index_sequence into calls to func
  template <class F, std::size_t... Is>
  XARA_UNROLL_INLINE constexpr void 
  repeat_impl(F func, std::index_sequence<Is...>) noexcept
  {
    (func(num<Is>{}), ...);
  }

  // Helper that unpacks the index_sequence into calls to func, offset by Start
  template <std::size_t Start, class F, std::size_t... Is>
  XARA_UNROLL_INLINE constexpr void 
  unroll_impl(F func, std::index_sequence<Is...>) noexcept
  {
      (func(num<Start + Is>{}), ...);
  }
}

#if 0 //  __cplusplus >= 202002L
template <std::size_t N, class F>
[[gnu::always_inline]] inline constexpr void
Repeat(F func) noexcept
{
  // Lambda with templated parameter pack (C++20 feature)
  ([]<std::size_t... Is>(F func, std::index_sequence<Is...>){
      (func(num<Is>{}), ...);
  })(func, std::make_index_sequence<N>{});
}
#else
template <std::size_t N, class F>
[[gnu::always_inline]] inline constexpr void
Repeat(F func) noexcept
{
  repeat_impl(func, std::make_index_sequence<N>{});
}
#endif 


// Unroll<3, 3>(func); // generates nothing
// Unroll<3, 4>(func); // generates 3
// Unroll<3, 7>(func); // generates 3,4,5,6
template <std::size_t Start, std::size_t Stop, class F>
inline constexpr void
Unroll(F func) noexcept
{
  static_assert(Stop >= Start, "Stop must be greater than or equal to Start");
  constexpr std::size_t N = Stop - Start;
  unroll_impl<Start>(func, std::make_index_sequence<N>{});
}

} // namespace OpenSees