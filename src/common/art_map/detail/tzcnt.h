#ifndef ART_DETAIL_TZCNT_HEADER_INCLUDED
#define ART_DETAIL_TZCNT_HEADER_INCLUDED

namespace art
{
namespace detail
{

// Specializations for the ctz builtin. Should we want to extend to other
// compilers than clang and g++, this would be the place to show one's preprocessor-fu.
template <unsigned N> struct tzcnt_selector;

template <> struct tzcnt_selector<sizeof(unsigned int)> {
    static int apply(unsigned int x) noexcept { return __builtin_ctz(x); }
};

template <> struct tzcnt_selector<sizeof(unsigned long)> {
    static int apply(unsigned long x) noexcept { return __builtin_ctzl(x); }
};

template <typename T> inline int tzcnt(T x) noexcept
{
    return tzcnt_selector<sizeof(T)>::apply(x);
}

} // namespace detail
} // namespace art

#endif // ART_DETAIL_TZCNT_HEADER_INCLUDED
