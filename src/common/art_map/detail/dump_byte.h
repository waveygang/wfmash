#ifndef ART_DETAIL_DUMP_BYTE_HEADER_INCLUDED
#define ART_DETAIL_DUMP_BYTE_HEADER_INCLUDED

#include <array>
#include <ostream>

namespace art
{
namespace detail
{

inline void dump_byte(std::ostream& os, std::uint8_t byte)
{
    char hex[4];
    std::snprintf(hex, 4, " %.2x", static_cast<unsigned>(byte));
    os.write(hex, 3);
}

// Prints any POD as bytes
template <typename T> inline void dump_bytes(std::ostream& os, T u)
{
    union {
        T value;
        std::array<std::uint8_t, sizeof(T)> bytes;
    } const val{u};
    for (std::uint8_t b : val.bytes)
        dump_byte(os, b);
}

} // namespace detail
} // namespace art

#endif // ART_DETAIL_DUMP_BYTE_HEADER_INCLUDED
