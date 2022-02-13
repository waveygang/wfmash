#ifndef ART_DETAIL_ART_NODE_BASE_HEADER_INCLUDED
#define ART_DETAIL_ART_NODE_BASE_HEADER_INCLUDED

#include "dump_byte.h"

#include <cassert>
#include <cstdint>

namespace art
{
namespace detail
{

// Common base for various node types and leaves. This common base ensures that
// common functionality could be achieved though the base pointer.
template <typename BitwiseKey> struct art_node_base {
    using bitwise_key = BitwiseKey;
    using key_size_type = typename BitwiseKey::size_type;

    explicit constexpr art_node_base(bitwise_key key) noexcept
        : header(key)
    {
    }

    [[nodiscard]] constexpr bitwise_key prefix() const noexcept { return header; }
    [[nodiscard]] constexpr std::uint8_t front() const noexcept { return header.front(); }

    // This is meaningful for internal nodes only
    [[nodiscard]] constexpr key_size_type prefix_length() const noexcept { return header.size(); }
    [[nodiscard]] constexpr key_size_type shared_prefix_length(bitwise_key key) const noexcept
    {
        return bitwise_key::shared_len(key, header, header.size());
    }
    [[nodiscard]] bitwise_key shared_prefix(bitwise_key key) const noexcept
    {
        return bitwise_key::partial_key(key, shared_prefix_length(key));
    }

    // Prefix manipulation routines
    constexpr void shift_right(key_size_type size) noexcept { header.shift_right_resize(size); }
    constexpr void shift_left(std::uint8_t key) noexcept { header.shift_left_resize(key); }
    constexpr void shift_left(bitwise_key key) noexcept { header.shift_left_resize(key); }

    // Dump bitwise key prefix
    static void dump(std::ostream& os, bitwise_key key, key_size_type len)
    {
        for (key_size_type i = 0; i != len; ++i)
            dump_byte(os, key[i]);
    }
    static void dump(std::ostream& os, bitwise_key key)
    {
        os << "key prefix len = " << static_cast<std::size_t>(key.size());
        if (key.size()) {
            os << ", key prefix =";
            dump(os, key, key.size());
        }
    }

private:
    BitwiseKey header;
};

} // namespace detail
} // namespace art

#endif // ART_DETAIL_ART_NODE_BASE_HEADER_INCLUDED
