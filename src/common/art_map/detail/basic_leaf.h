#ifndef ART_DETAIL_BASIC_LEAF_HEADER_INCLUDED
#define ART_DETAIL_BASIC_LEAF_HEADER_INCLUDED

#include "art_node_base.h"

namespace art
{
namespace detail
{

// Helper struct for leaf node-related data
template <typename BitwiseKey, typename T>
struct basic_leaf final : public art_node_base<BitwiseKey> {
    using value_type = T;
    using bitwise_key = BitwiseKey;
    using parent_type = art_node_base<BitwiseKey>;
    using key_type = typename BitwiseKey::key_type;

    template <typename... Args>
    constexpr basic_leaf(bitwise_key key,
                         Args&&... args) noexcept(std::is_nothrow_constructible<T>::value)
        : parent_type(key)
        , data(std::forward<Args>(args)...)
    {
    }

    // There is always a single element in this leaf
    [[nodiscard]] static constexpr std::size_t size() noexcept { return 1; }

    [[nodiscard]] T& value() noexcept { return data; }
    [[nodiscard]] const T& value() const noexcept { return data; }

    // Overwrite currently held value
    template <typename... Args>
    void push_front(Args&&... args) noexcept(
        std::is_nothrow_constructible<T>::value&& std::is_nothrow_move_assignable<T>::value)
    {
        T tmp(std::forward<Args>(args)...);
        data = std::move(tmp);
    }

    [[noreturn]] static void push_back(T&&) { throw std::runtime_error("basic_leaf: push_back"); }

    void dump(std::ostream& os) const
    {
        os << "key =";
        parent_type::dump(os, this->prefix(), this->prefix().max_size());
        os << ", value = " << value();
    }

private:
    T data;
};

// Specialization for integral constants that do not take any space
template <typename BitwiseKey, typename T, T V>
struct basic_leaf<BitwiseKey, std::integral_constant<T, V>> final
    : public art_node_base<BitwiseKey> {
    using value_type = std::integral_constant<T, V>;
    using bitwise_key = BitwiseKey;
    using parent_type = art_node_base<BitwiseKey>;
    using key_type = typename BitwiseKey::key_type;

    constexpr basic_leaf(bitwise_key key, value_type) noexcept
        : parent_type(key)
    {
    }

    // There is always a single element in this leaf
    [[nodiscard]] static constexpr std::size_t size() noexcept { return 1; }

    // Simply return a value
    [[nodiscard]] static constexpr value_type value() noexcept { return value_type(); }

    static constexpr void push_front(value_type) noexcept {}
    [[noreturn]] static void push_back(value_type)
    {
        throw std::runtime_error("basic_leaf: push_back");
    }

    void dump(std::ostream& os) const
    {
        os << "key =";
        parent_type::dump(os, this->prefix(), this->prefix().max_size());
    }
};

} // namespace detail
} // namespace art

#endif // ART_DETAIL_BASIC_LEAF_HEADER_INCLUDED
