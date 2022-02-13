#ifndef ART_DETAIL_ART_DELETER_HEADER_INCLUDED
#define ART_DETAIL_ART_DELETER_HEADER_INCLUDED

#include <memory>

namespace art
{
namespace detail
{

template <typename Node, typename Db> class node_deleter final
{
public:
    explicit node_deleter(Db& db) noexcept
        : db_{&db}
    {
    }
    void operator()(Node* node) const noexcept { db_->deallocate(node); }

private:
    Db* db_;
};

template <typename Node, typename Db>
using unique_node_ptr = std::unique_ptr<Node, node_deleter<Node, Db>>;

template <typename Node, typename Db>
inline unique_node_ptr<Node, Db> make_unique_node_ptr(Node* node, Db& db) noexcept
{
    return unique_node_ptr<Node, Db>(node, node_deleter<Node, Db>(db));
}

} // namespace detail
} // namespace art

#endif // ART_DETAIL_ART_DELETER_HEADER_INCLUDED
