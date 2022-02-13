#ifndef ART_DETAIL_NODE_TYPE_HEADER_INCLUDED
#define ART_DETAIL_NODE_TYPE_HEADER_INCLUDED

#include <cstdint>

namespace art
{
namespace detail
{

enum class node_type : std::uint8_t { LEAF, I4, I16, I64, I256 };

} // namespace detail
} // namespace art

#endif // ART_DETAIL_NODE_TYPE_HEADER_INCLUDED
