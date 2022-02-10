#ifndef ART_DETAIL_ART_CONTAINER_IMPL_HEADER_INCLUDED
#define ART_DETAIL_ART_CONTAINER_IMPL_HEADER_INCLUDED

#include "art_container.h"
#include "art_deleter.h"
#include "art_nodes.h"

#include <boost/config.hpp>

namespace art
{
namespace detail
{

template <typename BitwiseKey>
inline constexpr void shift_right(std::pair<BitwiseKey, typename BitwiseKey::size_type>& lhs,
                                  typename BitwiseKey::size_type size) noexcept
{
    assert(lhs.second >= size);
    lhs.first.shift_right(size);
    lhs.second -= size;
}

template <typename P>
inline typename db<P>::const_iterator db<P>::internal_locate(bitwise_key_prefix& key) const noexcept
{
    // By ensuring that the leaf tag is always 0, we can simplify inode filtering
    // because nullptr pointers would serendipitously have 0 tags too
    static_assert(static_cast<unsigned>(node_type::LEAF) == 0, "Leaf tag must be 0");

    const_iterator pos(tree.root);

    while (pos.tag() != node_type::LEAF) {
        const inode* nb = pos.inode();
        const key_size_type prefix_length = nb->prefix_length();
        if (prefix_length &&
            (key.second <= prefix_length || nb->shared_prefix_length(key.first) < prefix_length))
            break;

        const auto child = inode::find_child(pos.node(), key.first[prefix_length]);
        if (!child.node())
            break;
        pos = child;

        // Consume the explored prefix + 1 byte used during child lookup
        shift_right(key, prefix_length + 1);
    }

    return pos;
}

template <typename BitwiseKey>
[[nodiscard]] inline constexpr auto make_key_prefix(BitwiseKey key) noexcept
{
    return std::make_pair(key, key.max_size());
}

template <typename P>
inline typename db<P>::const_iterator db<P>::internal_find(fast_key_type key) const noexcept
{
    const bitwise_key bitk(key);
    auto pos = internal_locate(make_key_prefix(bitk));
    return pos.match(bitk) ? pos : end();
}

template <typename P> inline bool db<P>::contains(fast_key_type key) const noexcept
{
    const bitwise_key bitk(key);
    return internal_locate(make_key_prefix(bitk)).match(bitk);
}

template <typename P>
inline typename db<P>::const_iterator db<P>::forward_step(const_iterator pos) const noexcept
{
    return pos.parent() != nullptr ? pos.forward_step() : end();
}

template <typename P>
template <typename Filter>
inline typename db<P>::const_iterator db<P>::internal_bound(bitwise_key original_key,
                                                            Filter filter) const noexcept
{
    // By ensuring that the leaf tag is always 0, we can simplify inode filtering
    // because nullptr pointers would serendipitously have 0 tags too
    static_assert(static_cast<unsigned>(node_type::LEAF) == 0, "Leaf tag must be 0");

    if (BOOST_UNLIKELY(empty())) {
        return end();
    }

    const_iterator pos(tree.root);

    auto key = make_key_prefix(original_key);

    while (pos.tag() != node_type::LEAF) {
        const inode* nb = pos.inode();
        const key_size_type prefix_length = nb->prefix_length();
        if (prefix_length &&
            (key.second <= prefix_length || nb->shared_prefix_length(key.first) < prefix_length)) {
            return forward_step(pos);
        }

        const std::uint8_t key_byte = key.first[prefix_length];
        const auto p = inode::lower_bound(pos.node(), key_byte);
        if (!p.first.node()) {
            return forward_step(pos);
        }
        if (p.second > key_byte) {
            return p.first.tag() != node_type::LEAF ? inode::leftmost_leaf(p.first.node())
                                                    : p.first;
        }
        pos = p.first;

        // Consume the explored prefix + 1 byte used during child lookup
        shift_right(key, prefix_length + 1);
    }

    // We found some leaf. At this point we know that the prefix of that leaf
    // matched the given key prefix, but we don't know much about the concrete
    // value in that leaf. The value might well be smaller than the one we seek.
    return filter(pos.leaf()->prefix().unpack()) ? forward_step(pos) : pos;
}

template <typename P>
inline typename db<P>::const_iterator db<P>::internal_lower_bound(fast_key_type key) const noexcept
{
    return internal_bound(bitwise_key(key),
                          [key, comp{key_comp()}](fast_key_type leaf) { return comp(leaf, key); });
}

template <typename P>
inline typename db<P>::const_iterator db<P>::internal_upper_bound(fast_key_type key) const noexcept
{
    return internal_bound(bitwise_key(key),
                          [key, comp{key_comp()}](fast_key_type leaf) { return !comp(key, leaf); });
}

template <typename P>
template <typename Node, typename... Args>
inline unique_node_ptr<Node, db<P>> db<P>::make_node_ptr(Args&&... args)
{
    using unique_ptr_t = unique_node_ptr<Node, db<P>>;

    using node_allocator_type = typename db_allocator_traits::template rebind_alloc<Node>;
    using node_allocator_traits = std::allocator_traits<node_allocator_type>;

    node_allocator_type alloc(allocator());
    unique_ptr_t node_ptr(node_allocator_traits::allocate(alloc, 1),
                          node_deleter<Node, self_t>(*this));

    node_allocator_traits::construct(alloc, node_ptr.get(), std::forward<Args>(args)...);

    ++count<Node>();

    return node_ptr;
}

template <typename P>
template <typename Node>
inline void db<P>::deallocate_node(Node* node) noexcept
{
    using node_allocator_type = typename db_allocator_traits::template rebind_alloc<Node>;
    using node_allocator_traits = std::allocator_traits<node_allocator_type>;

    assert(count<Node>() != 0);

    node_allocator_type alloc(allocator());
    node_allocator_traits::destroy(alloc, node);
    node_allocator_traits::deallocate(alloc, node, 1);

    --count<Node>();
}

template <typename P>
inline void db<P>::deallocate(leaf_type* leaf) noexcept(
    std::is_nothrow_destructible<mapped_type>::value)
{
    deallocate_node(leaf);
}

template <typename P>
template <typename Node>
inline void db<P>::deallocate_subtree(Node* node) noexcept(
    std::is_nothrow_destructible<mapped_type>::value)
{
    node->delete_subtree(*this);
    deallocate_node(node);
}

// Deallocate a tagged pointer
template <typename P>
inline void db<P>::deallocate(node_ptr node) noexcept(
    std::is_nothrow_destructible<mapped_type>::value)
{
    switch (node.tag()) {
    case node_type::LEAF:
        deallocate(static_cast<leaf_type*>(node.get()));
        break;
    case node_type::I4:
        deallocate_subtree(static_cast<inode_4*>(node.get()));
        break;
    case node_type::I16:
        deallocate_subtree(static_cast<inode_16*>(node.get()));
        break;
    case node_type::I64:
        deallocate_subtree(static_cast<inode_64*>(node.get()));
        break;
    case node_type::I256:
        deallocate_subtree(static_cast<inode_256*>(node.get()));
        break;
    default:
        ART_DETAIL_CANNOT_HAPPEN();
    }
}

template <typename P> inline void db<P>::assign_to_parent(iterator hint, node_ptr child) noexcept
{
    if (BOOST_LIKELY(hint.parent() != nullptr)) {
        const auto parent = hint.parent();
        inode::dispatch_inode(parent,
                              [hint, child](auto& inode) { return inode.replace(hint, child); });
    } else {
        tree.root = child;
    }
}

template <typename P>
template <typename INode, typename NodePtr, typename SizeType>
inline typename db<P>::iterator db<P>::create_inode(iterator hint, bitwise_key prefix, NodePtr pdst,
                                                    leaf_unique_ptr leaf, SizeType key)
{
    auto inode = make_node_ptr<INode>(prefix);
    const iterator leaf_iter = inode->populate(std::move(pdst), std::move(leaf), key);
    assign_to_parent(hint, leaf_iter.parent());
    inode.release(); // All went well, we can release the pointer
    return leaf_iter;
}

template <typename P>
template <typename Source>
inline typename db<P>::iterator db<P>::grow_node(iterator hint, leaf_unique_ptr leaf,
                                                 std::uint8_t key_byte)
{
    using larger_inode = typename Source::larger_inode_type;

    assert(hint.tag() == Source::static_type());

    auto dst = static_cast<Source*>(hint.inode());
    if (BOOST_LIKELY(!dst->is_full())) {
        return dst->add(std::move(leaf), key_byte);
    } else {
        // Destination node is full, needs to grow
        return create_inode<larger_inode>(hint, dst->prefix(), make_unique_node_ptr(dst, *this),
                                          std::move(leaf), key_byte);
    }
}

template <typename P>
template <typename... Args>
inline typename db<P>::iterator db<P>::internal_emplace(iterator hint, bitwise_key bitk,
                                                        const bitwise_key_prefix& key,
                                                        Args&&... args)
{
    assert(bitk.max_size() >= key.second);

    // Preemptively create a leaf. This also ensures strong exception safety
    auto leaf_ptr = make_node_ptr<leaf_type>(bitk, std::forward<Args>(args)...);

    if (BOOST_UNLIKELY(empty())) {
        assert(!hint.node());
        tree.root = node_ptr(leaf_ptr.release(), node_type::LEAF);
        return iterator(tree.root);
    }

    assert(hint.node());
    assert(key.second != 0);

    const node_type dst_tag = hint.tag();

    if (dst_tag == node_type::LEAF) {
        leaf_type* const pdst = hint.leaf();

        bitwise_key prefix = pdst->prefix();

        // Can only happen in multivalued container case
        if (BOOST_UNLIKELY(prefix == bitk)) {
            pdst->push_back(std::move(leaf_ptr->value()));
            return iterator(hint);
        }

        const key_size_type depth = bitk.max_size() - key.second;

        // Put the 2 leaves under the single inode_4
        prefix.shift_right(depth);
        const key_size_type len = bitwise_key::shared_len(prefix, key.first, key.second - 1);

        return create_inode<inode_4>(hint, bitwise_key::partial_key(prefix, len), pdst,
                                     std::move(leaf_ptr), depth);
    }

    // Some other node, not a leaf
    inode* const pdst = hint.inode();
    const key_size_type prefix_len = pdst->prefix_length();

    // If the prefix is not fully shared, split the prefix
    {
        const bitwise_key shared_prefix = pdst->shared_prefix(key.first);
        if (shared_prefix.size() < prefix_len) {
            return create_inode<inode_4>(hint, shared_prefix, hint.node(), std::move(leaf_ptr),
                                         key.first[shared_prefix.size()]);
        }
    }

    assert(key.second > prefix_len);

    // Key byte for this leaf. Pass it on to the nodes
    const std::uint8_t key_byte = key.first[prefix_len];

    // Add the newly created leaf to the proper node.
    // If that node is full, then a new larger node will be created,
    // and the leaf will be added there.
    switch (dst_tag) {
    case node_type::I4:
        return grow_node<inode_4>(hint, std::move(leaf_ptr), key_byte);
    case node_type::I16:
        return grow_node<inode_16>(hint, std::move(leaf_ptr), key_byte);
    case node_type::I64:
        return grow_node<inode_64>(hint, std::move(leaf_ptr), key_byte);
    case node_type::I256:
        return static_cast<inode_256*>(pdst)->add(std::move(leaf_ptr), key_byte);
    default:
        ART_DETAIL_CANNOT_HAPPEN();
    }
}

template <typename P>
template <typename... Args>
inline typename db<P>::iterator db<P>::emplace_key_args(std::true_type, bitwise_key key,
                                                        Args&&... args)
{
    // Emplacement support for multiset/multimap
    auto prefix = make_key_prefix(key);
    auto pos = iterator(internal_locate(prefix));
    return internal_emplace(pos, key, prefix, std::forward<Args>(args)...);
}

// Inserts a value into the tree only if it does not already exist. The
// boolean return value indicates whether insertion succeeded or failed.
template <typename P>
template <typename... Args>
inline std::pair<typename db<P>::iterator, bool> db<P>::emplace_key_args(std::false_type,
                                                                         bitwise_key key,
                                                                         Args&&... args)
{
    // Emplacement support for set/map
    auto prefix = make_key_prefix(key);
    auto pos = iterator(internal_locate(prefix));

    const bool should_insert = !pos.match(key);
    return std::make_pair(
        should_insert ? internal_emplace(pos, key, prefix, std::forward<Args>(args)...) : pos,
        should_insert);
}

template <typename Iterator>
[[nodiscard]] inline constexpr Iterator get_iterator(Iterator it) noexcept
{
    return it;
}

template <typename Iterator>
[[nodiscard]] inline constexpr Iterator get_iterator(std::pair<Iterator, bool> p) noexcept
{
    return p.first;
}

template <typename P>
template <typename... Args>
inline typename db<P>::iterator db<P>::emplace_hint_key_args(iterator hint, fast_key_type key,
                                                             Args&&... args)
{
    const bitwise_key bitk(key);
    // Fast path for insert(lower_bound())
    if (hint.match(bitk)) {
        hint.leaf()->push_front(std::forward<Args>(args)...);
        return hint;
    }
    return get_iterator(emplace_key_args(multi_container(), bitk, std::forward<Args>(args)...));
}

template <typename P>
template <typename Source>
inline typename db<P>::iterator db<P>::shrink_node(iterator pos)
{
    const auto node = pos.parent();
    assert(node.tag() == Source::static_type());

    auto src = static_cast<Source*>(node.get());
    assert(pos.parent() == src->tagged_self());

    if (BOOST_LIKELY(!src->is_min_size())) {
        src->remove(pos.index());
    } else {
        // If allocation fails here, the original node will be left untouched,
        // which gives the shrinking operation the strong exception safety guarantee
        auto smaller = make_tagged_ptr(make_smaller_node(*src, pos));
        assign_to_parent(src->self_iterator(), smaller);
        // All went well, we can deallocate the original node
        deallocate(src);
    }
    pos.position -= 1;
    return pos;
}

template <typename P> inline typename db<P>::iterator db<P>::internal_erase(iterator pos)
{
    assert(pos.is_leaf());

    iterator after_erase;

    // Also remove a leaf from its parent. In case of the no parent case,
    // we'll simply delete the root leaf
    const auto leaf_parent = pos.parent();

    if (BOOST_UNLIKELY(leaf_parent == nullptr)) {
        assert(pos.node() == tree.root);
        tree.root = nullptr;
        after_erase = iterator{};
    } else {
        // Remove the value from to the parent node.
        // If that node is already of minimum size, then a new smaller
        // node will be created without the value to remove.
        switch (leaf_parent.tag()) {
        case node_type::I4:
            after_erase = shrink_node<inode_4>(pos);
            break;
        case node_type::I16:
            after_erase = shrink_node<inode_16>(pos);
            break;
        case node_type::I64:
            after_erase = shrink_node<inode_64>(pos);
            break;
        case node_type::I256:
            after_erase = shrink_node<inode_256>(pos);
            break;
        default:
            ART_DETAIL_CANNOT_HAPPEN();
        }
    }

    // All went well, we can deallocate the leaf
    deallocate(pos.leaf());
    return after_erase;
}

template <typename P>
inline typename db<P>::size_type db<P>::count(fast_key_type key) const noexcept
{
    const bitwise_key bitk(key);
    auto pos = internal_locate(make_key_prefix(bitk));
    return pos.match(bitk) ? pos.leaf()->size() : static_cast<size_type>(0);
}

template <typename P> inline typename db<P>::size_type db<P>::erase(fast_key_type key)
{
    const bitwise_key bitk(key);
    auto pos = internal_locate(make_key_prefix(bitk));

    size_type removed = 0;
    if (pos.match(bitk)) {
        removed = pos.leaf()->size();
        internal_erase(iterator(pos));
    }
    return removed;
}

template <typename P> inline void db<P>::swap(self_t& other) noexcept
{
    // Swap the tree
    std::swap(allocator(), other.allocator());
    std::swap(tree.root, other.tree.root);

    // Swap the stats
    std::swap(count_, other.count_);
}

template <typename P> inline void db<P>::clear()
{
    if (!empty()) {
        deallocate(tree.root);
        assert(leaf_count() == 0);
        tree.root = nullptr;
    }
}

template <typename P> inline void db<P>::dump(std::ostream& os) const
{
    os << "leaf size = " << sizeof(leaf_type) << ", inode size = " << sizeof(inode)
       << ", current memory use = " << current_memory_use() << '\n';
    inode::dump(os, tree.root);
}

} // namespace detail
} // namespace art

#endif // ART_DETAIL_ART_CONTAINER_IMPL_HEADER_INCLUDED
