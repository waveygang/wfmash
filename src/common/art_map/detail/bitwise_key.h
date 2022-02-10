#ifndef ART_DETAIL_BITWISE_KEY_HEADER_INCLUDED
#define ART_DETAIL_BITWISE_KEY_HEADER_INCLUDED

#include "tzcnt.h"

#include <boost/endian/conversion.hpp>
#include <boost/integer.hpp>

#include <cassert>
#include <climits>
#include <functional>
#include <type_traits>

namespace art
{
namespace detail
{

namespace comparison_ops
{

struct less_tag {
};

struct greater_tag {
};

// Internal ART key in binary-comparable format
template <typename T> struct compared_arg;

// Specializations for particular orderings
template <typename T> struct compared_arg<std::less<T>> {
    using key_type = T;
    using order = less_tag;
};

template <typename T> struct compared_arg<std::greater<T>> {
    using key_type = T;
    using order = greater_tag;
};

} // namespace comparison_ops

template <typename Key, typename Compare> struct bitwise_compare;

// Specializations for particular orderings
template <typename Key> struct bitwise_compare<Key, comparison_ops::less_tag> {
    using key_type = Key;
    using bitkey_type = Key;

    [[nodiscard]] inline static constexpr Key byte_swap(Key k) noexcept
    {
        // Flip bytes to Big Endian order (no-op on Big Endian architectures)
        return boost::endian::native_to_big(k);
    }
    [[nodiscard]] inline static constexpr Key unpack(Key k) noexcept
    {
        return boost::endian::big_to_native(k);
    }
};

template <typename Key> struct bitwise_compare<Key, comparison_ops::greater_tag> {
    using key_type = Key;
    using bitkey_type = Key;

    [[nodiscard]] inline static constexpr Key byte_swap(Key k) noexcept
    {
        return boost::endian::native_to_big(~k);
    }
    [[nodiscard]] inline static constexpr Key unpack(Key k) noexcept
    {
        return ~boost::endian::big_to_native(k);
    }
};

template <typename Ptr, typename Order> struct ptr_bitwise_compare {
    static_assert(std::is_pointer<Ptr>::value, "Unsupported pointer type");

    using key_type = Ptr;
    using bitkey_type = std::uintptr_t;
    using compare_t = bitwise_compare<bitkey_type, Order>;

    [[nodiscard]] inline static constexpr bitkey_type byte_swap(Ptr k) noexcept
    {
        return compare_t::byte_swap(reinterpret_cast<bitkey_type>(k));
    }
    [[nodiscard]] inline static constexpr Ptr unpack(bitkey_type k) noexcept
    {
        return reinterpret_cast<Ptr>(compare_t::unpack(k));
    }
};

template <typename Int, typename Order> struct int_bitwise_compare {
    static_assert(std::is_signed<Int>::value && !std::is_floating_point<Int>::value,
                  "Unsupported signed integer type");

    using key_type = Int;
    using bitkey_type = typename boost::uint_t<sizeof(Int) * CHAR_BIT>::fast;
    using compare_t = bitwise_compare<bitkey_type, Order>;

    static constexpr std::size_t num_digits = std::numeric_limits<bitkey_type>::digits;
    static constexpr bitkey_type sign_bit = static_cast<bitkey_type>(1) << (num_digits - 1);

    [[nodiscard]] inline static constexpr bitkey_type byte_swap(Int k) noexcept
    {
        return compare_t::byte_swap(static_cast<bitkey_type>(k) ^ sign_bit);
    }
    [[nodiscard]] inline static constexpr Int unpack(bitkey_type k) noexcept
    {
        return static_cast<Int>(compare_t::unpack(k) ^ sign_bit);
    }
};

template <typename Compare> struct unsigned_integral_bitwise_key {

    using bitkey_type = typename Compare::bitkey_type;
    static_assert(std::is_unsigned<bitkey_type>::value && std::is_integral<bitkey_type>::value,
                  "Unsupported unsigned integral key type");

    using key_type = typename Compare::key_type;
    static_assert(sizeof(bitkey_type) == sizeof(key_type), "Invalid key size");

    // std::uint8_t would be sufficient here, but it leads to poorer
    // performance in benchmarks.
    using size_type = unsigned int;
    static constexpr size_type num_bytes = sizeof(key_type);

    [[nodiscard]] static constexpr size_type max_size() noexcept { return num_bytes; }

    constexpr unsigned_integral_bitwise_key() noexcept = default;
    explicit constexpr unsigned_integral_bitwise_key(key_type k) noexcept
        : key{Compare::byte_swap(k)}
    {
    }

    [[nodiscard]] size_type size() const noexcept { return key.bytes[num_bytes - 1]; }

    // Internal nodes use partial keys
    [[nodiscard]] std::uint8_t operator[](size_type index) const noexcept
    {
        assert(index < num_bytes);
        return key.bytes[index];
    }

    bool operator==(unsigned_integral_bitwise_key rhs) const noexcept
    {
        return key.bits == rhs.key.bits;
    }

    [[nodiscard]] constexpr std::uint8_t front() const noexcept { return key.bytes[0]; }

    constexpr void shift_right(size_type nbytes) noexcept { key.bits >>= (nbytes * CHAR_BIT); }
    constexpr void shift_right_resize(size_type nbytes) noexcept
    {
        assert(nbytes <= size());
        const size_type new_size = size() - nbytes;
        shift_right(nbytes);
        put_size(new_size);
    }
    constexpr void shift_left_resize(std::uint8_t value) noexcept
    {
        const size_type new_size = size() + 1;
        assert(new_size < max_size());
        shift_left(1);
        key.bytes[0] = value;
        put_size(new_size);
    }
    constexpr void shift_left_resize(unsigned_integral_bitwise_key value) noexcept
    {
        const size_type len = value.size();
        const size_type new_size = size() + len;
        assert(new_size < max_size());
        shift_left(len);
        key.bits |= value.key.bits;
        put_size(new_size);
    }

    [[nodiscard]] static size_type shared_len(unsigned_integral_bitwise_key k1,
                                              unsigned_integral_bitwise_key k2,
                                              size_type clamp_byte_pos) noexcept
    {
        assert(clamp_byte_pos < num_bytes);

        const bitkey_type diff = k1.key.bits ^ k2.key.bits;
        const bitkey_type clamped = diff | himask(clamp_byte_pos);
        return tzcnt(clamped) >> 3U;
    }

    [[nodiscard]] static unsigned_integral_bitwise_key partial_key(unsigned_integral_bitwise_key k,
                                                                   size_type cut_len) noexcept
    {
        unsigned_integral_bitwise_key p(k.key.bits & (himask(cut_len) - 1), std::false_type());
        p.put_size(cut_len);
        return p;
    }

    [[nodiscard]] key_type unpack() const noexcept { return Compare::unpack(key.bits); }

private:
    // Non-byte-swapping constructor
    constexpr unsigned_integral_bitwise_key(bitkey_type k, std::false_type) noexcept
        : key{k}
    {
    }

    static constexpr bitkey_type himask(size_type len) noexcept
    {
        return static_cast<bitkey_type>(1) << (len * CHAR_BIT);
    }

    constexpr void shift_left(size_type nbytes) noexcept { key.bits <<= (nbytes * CHAR_BIT); }
    constexpr void put_size(size_type len) noexcept { key.bytes[num_bytes - 1] = len; }

    // Define a SIMD byte vector. This has an advantage over, say,
    // std::array<std::uint8_t, num_bytes> that it tells the compiler in no uncertain terms that
    // this thing in a SIMD vector and should be treated as such. Let the compiler do the heavy
    // lifting how to optimize this beast.
    typedef std::uint8_t byte_vec __attribute__((vector_size(sizeof(std::uint8_t) * num_bytes)));

    union {
        bitkey_type bits;
        byte_vec bytes;
    } key;
};

// Support for unsigned keys
template <typename Key, typename Order> struct unsigned_bitwise_key {
    using type = unsigned_integral_bitwise_key<bitwise_compare<Key, Order>>;
};

template <typename Key, typename Order>
using unsigned_bitwise_key_t = typename unsigned_bitwise_key<Key, Order>::type;

// Supports for non member-function pointers
template <typename Key, typename Order> struct pointer_bitwise_key {
    using type = unsigned_integral_bitwise_key<ptr_bitwise_compare<Key, Order>>;
};

template <typename Key, typename Order>
using pointer_bitwise_key_t = typename pointer_bitwise_key<Key, Order>::type;

// Supports for signed integers of various sizes. Floats and doubles are
// explicitly not supported
template <typename Key, typename Order> struct signed_bitwise_key {
    using type = unsigned_integral_bitwise_key<int_bitwise_compare<Key, Order>>;
};

template <typename Key, typename Order>
using signed_bitwise_key_t = typename signed_bitwise_key<Key, Order>::type;

template <typename Key, typename Compare> struct bitwise_key_compare {
    using predicate = comparison_ops::compared_arg<Compare>;

    using arg_t = typename predicate::key_type;
    static_assert(std::is_same<Key, arg_t>::value, "Incompatible comparison predicate");

    using order = typename predicate::order;
    using type = std::conditional_t<
        std::is_signed<arg_t>::value, signed_bitwise_key_t<arg_t, order>,
        std::conditional_t<std::is_pointer<arg_t>::value, pointer_bitwise_key_t<arg_t, order>,
                           unsigned_bitwise_key_t<arg_t, order>>>;
};

template <typename Key, typename Compare>
using bitwise_key_t = typename bitwise_key_compare<Key, Compare>::type;

} // namespace detail
} // namespace art

#endif // ART_DETAIL_BITWISE_KEY_HEADER_INCLUDED
