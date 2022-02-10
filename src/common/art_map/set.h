#ifndef ART_SET_HEADER_INCLUDED
#define ART_SET_HEADER_INCLUDED

#include "detail/art_container.h"
#include "detail/container_traits.h"

namespace art
{
namespace detail
{

// Set traits
template <typename Key, typename Compare, typename Alloc, typename MultiSet>
struct set_traits : public container_traits<Key, std::false_type, Key, Compare, Alloc, MultiSet> {
    using mapped_type = std::false_type;

    using fast_key_type = fast_const_argument_t<Key>;
    using reference = fast_key_type;
    using const_reference = fast_key_type;

    [[nodiscard]] static fast_key_type key(fast_key_type x) noexcept { return x; }

    template <typename T> [[nodiscard]] static constexpr mapped_type value(T&&) noexcept
    {
        return mapped_type();
    }

    [[nodiscard]] static fast_key_type value_ref(Key key, mapped_type) noexcept { return key; }
};

} // namespace detail

// Parametrization for unique element set
template <typename Key, typename Compare = std::less<Key>, typename Alloc = std::allocator<Key>>
using set = detail::db<detail::set_traits<Key, Compare, Alloc, std::false_type>>;

// Parametrization for multiset
template <typename Key, typename Compare = std::less<Key>, typename Alloc = std::allocator<Key>>
using multiset = detail::db<detail::set_traits<Key, Compare, Alloc, std::true_type>>;

} // namespace art

#endif // ART_SET_HEADER_INCLUDED
