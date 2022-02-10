#ifndef ART_DETAIL_TAGGED_PTR_H
#define ART_DETAIL_TAGGED_PTR_H

#include <boost/integer/integer_mask.hpp>
#include <boost/integer/static_log2.hpp>

#include <cstdint>
#include <stdexcept>

namespace art
{
namespace detail
{
namespace tagged
{

// Direct policy stores the tag within the managed pointer itself. Custom policies can be provided
template <typename T, unsigned int Bits = boost::static_log2<alignof(T) * CHAR_BIT>::value,
          typename TagType = typename boost::uint_t<Bits>::fast>
struct direct final {
    static_assert(Bits != 0 && Bits <= boost::static_log2<alignof(T) * CHAR_BIT>::value,
                  "Invalid tag bits for type T");

    using tag_type = TagType;

    static constexpr std::uintptr_t tag_mask = boost::low_bits_mask_t<Bits>::sig_bits_fast;
    static constexpr std::uintptr_t ptr_mask = ~tag_mask;

    [[nodiscard]] static std::uintptr_t uintptr(const T* p) noexcept
    {
        return reinterpret_cast<std::uintptr_t>(p);
    }

    [[nodiscard]] static std::uintptr_t ptr_bits(const T* p) noexcept
    {
        return uintptr(p) & ptr_mask;
    }

    [[nodiscard]] static T* extract(T* p) noexcept { return reinterpret_cast<T*>(ptr_bits(p)); }

    [[nodiscard]] static tag_type tag(const T* p) noexcept
    {
        return static_cast<tag_type>(uintptr(p) & tag_mask);
    }

    [[nodiscard]] static constexpr T* tag(const T* p, tag_type value) noexcept
    {
        // Take care that the tag value does not step on pointer bits
        assert(!(ptr_mask & static_cast<std::uintptr_t>(value)));
        return reinterpret_cast<T*>(ptr_bits(p) | static_cast<std::uintptr_t>(value));
    }
};

} // namespace tagged

// Tagged pointer is a pointer that allows to access lower memory bits in the pointer value itself.
// These lower memory bits are unused due to alignment, so there we can put some information tied to
// the pointer that would otherwise require an external container.
template <typename T, typename Policy = tagged::direct<T>> class tagged_ptr final
{
public:
    using tag_type = typename Policy::tag_type;
    using pointer = T*;
    using const_pointer = std::add_const_t<pointer>;
    using reference = T&;

    tagged_ptr() noexcept = default;
    tagged_ptr(pointer ptr, tag_type init) noexcept
        : ptr_(Policy::tag(ptr, init))
    {
    }

    [[nodiscard]] pointer get() const noexcept { return Policy::extract(ptr_); }

    explicit operator bool() const noexcept { return ptr_ != nullptr; }
    bool operator==(std::nullptr_t) const noexcept { return ptr_ == nullptr; }
    bool operator!=(std::nullptr_t) const noexcept { return ptr_ != nullptr; }
    bool operator==(tagged_ptr rhs) const noexcept { return ptr_ == rhs.ptr_; }
    bool operator!=(tagged_ptr rhs) const noexcept { return ptr_ != rhs.ptr_; }
    tagged_ptr& operator=(std::nullptr_t) noexcept
    {
        ptr_ = nullptr;
        return *this;
    }
    tagged_ptr& operator=(pointer p) noexcept // Zero tag assignment
    {
        ptr_ = p;
        return *this;
    }

    reference operator*() const noexcept { return *get(); }
    pointer operator->() noexcept { return get(); }
    const_pointer operator->() const noexcept { return get(); }

    [[nodiscard]] tag_type tag() const noexcept { return Policy::tag(ptr_); }
    void tag(tag_type value) { ptr_ = Policy::tag(ptr_, value); }

    [[nodiscard]] static tag_type tag(const_pointer p) noexcept { return Policy::tag(p); }

private:
    pointer ptr_;
};

} // namespace detail
} // namespace art

#endif // ART_DETAIL_TAGGED_PTR_H
