#ifndef ART_DETAIL_ART_CONTAINER_HEADER_INCLUDED
#define ART_DETAIL_ART_CONTAINER_HEADER_INCLUDED

#include "art_deleter.h"
#include "tree_iterator.h"

#include <cassert>

namespace art
{
namespace detail
{
// Forward declarations for node types
template <typename Db> class basic_inode_impl;
template <typename Db> class basic_inode_4;
template <typename Db> class basic_inode_16;
template <typename Db> class basic_inode_64;
template <typename Db> class basic_inode_256;

template <typename Traits> class db
{
public:
    using bitwise_key = typename Traits::bitwise_key;

private:
    using multi_container = typename Traits::multi_container;
    using node_base = typename Traits::node_base;
    using node_ptr = typename Traits::node_ptr;
    using fast_key_type = typename Traits::fast_key_type;
    using leaf_type = typename Traits::leaf_type;

    using self_t = db<Traits>;
    using leaf_unique_ptr = unique_node_ptr<leaf_type, self_t>;

    using inode = basic_inode_impl<self_t>;
    using inode_4 = basic_inode_4<self_t>;
    using inode_16 = basic_inode_16<self_t>;
    using inode_64 = basic_inode_64<self_t>;
    using inode_256 = basic_inode_256<self_t>;

    // We will be friends only with the nodes that have the same policy
    friend inode;
    friend inode_4;
    friend inode_16;
    friend inode_64;
    friend inode_256;

    // Make deleters friends
    friend node_deleter<inode_4, self_t>;
    friend node_deleter<inode_16, self_t>;
    friend node_deleter<inode_64, self_t>;
    friend node_deleter<inode_256, self_t>;
    friend node_deleter<node_base, self_t>;
    friend node_deleter<leaf_type, self_t>;

public:
    using key_type = typename Traits::key_type;
    using mapped_type = typename Traits::mapped_type;
    using value_type = typename Traits::value_type;
    using pointer = typename Traits::pointer;
    using const_pointer = typename Traits::const_pointer;
    using reference = typename Traits::reference;
    using const_reference = typename Traits::const_reference;
    using size_type = typename Traits::size_type;
    using difference_type = typename Traits::difference_type;
    using key_compare = typename Traits::key_compare;
    using allocator_type = typename Traits::allocator_type;
    using iterator = tree_iterator<Traits, inode>;
    using const_iterator = tree_iterator<Traits, const inode>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

public:
    // Creation and destruction as defined in C++14 for std::map
    constexpr db() = default;

    explicit db(const allocator_type& alloc)
        : tree(alloc)
    {
    }

    template <class InputIt>
    db(InputIt first, InputIt last, const allocator_type& alloc)
        : tree(alloc)
    {
        insert(first, last);
    }

    db(const db& other);
    db(const db& other, const allocator_type& alloc);
    db(db&& other)
        : tree(std::move(other.tree))
        , count_(std::move(other.count_))
    {
        other.tree.root = nullptr;
    }
    db(db&& other, const allocator_type& alloc)
        : tree(alloc)
        , count_(std::move(other.count_))
    {
        tree.root = other.tree.root;
        other.tree.root = nullptr;
    }
    db(std::initializer_list<value_type> init, const allocator_type& alloc)
        : tree(alloc)
    {
        insert(init.begin(), init.end());
    }

    ~db() noexcept(std::is_nothrow_destructible<mapped_type>::value)
    {
        if (!empty())
            deallocate(tree.root);
    }

    // Allocator routines
    allocator_type get_allocator() const noexcept
    {
        return *static_cast<const allocator_type*>(&tree);
    }

    // We are not using the given comparator, so the function here is a dummy one
    key_compare key_comp() const { return key_compare(); }

    // Iterator routines.
    iterator begin() noexcept { return iterator(cbegin()); };
    const_iterator cbegin() const noexcept { return inode::leftmost_leaf(tree.root); }
    const_iterator begin() const noexcept { return cbegin(); }

    iterator end() noexcept { return iterator(tree.root, iterator::is_leaf(tree.root)); }
    const_iterator cend() const noexcept
    {
        return const_iterator(tree.root, iterator::is_leaf(tree.root));
    }
    const_iterator end() const noexcept { return cend(); }

    reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const noexcept { return const_reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const noexcept { return const_reverse_iterator(begin()); }

    size_type size() const noexcept { return leaf_count(); }
    bool empty() const noexcept { return tree.root == nullptr; }

    // Lookup routines
    iterator find(fast_key_type key) noexcept { return iterator(internal_find(key)); }
    const_iterator find(fast_key_type key) const noexcept { return internal_find(key); }

    iterator lower_bound(fast_key_type key) noexcept { return iterator(internal_lower_bound(key)); }
    const_iterator lower_bound(fast_key_type key) const noexcept
    {
        return internal_lower_bound(key);
    }

    iterator upper_bound(fast_key_type key) noexcept { return iterator(internal_upper_bound(key)); }
    const_iterator upper_bound(fast_key_type key) const noexcept
    {
        return internal_upper_bound(key);
    }

    // Since C++20
    bool contains(fast_key_type key) const noexcept;

    std::pair<iterator, iterator> equal_range(fast_key_type key)
    {
        return std::make_pair(lower_bound(key), upper_bound(key));
    }
    std::pair<const_iterator, const_iterator> equal_range(fast_key_type key) const
    {
        return std::make_pair(lower_bound(key), upper_bound(key));
    }

    // Insertion routines
    auto insert(const value_type& value)
    {
        return emplace_key_args(Traits::key(value), Traits::value(value));
    }
    auto insert(value_type&& value)
    {
        return emplace_key_args(Traits::key(value), Traits::value(std::move(value)));
    }
    template <class Arg> auto insert(Arg&& value)
    {
        return emplace_key_args(Traits::key(value), Traits::value(std::move(value)));
    }

    // Insert with hint. Check to see if the value should be placed immediately
    // before position in the tree. If it does, then the insertion will take
    // amortized constant time. If not, the insertion will take amortized
    // logarithmic time as if a call to insert(value) were made.
    iterator insert(iterator hint, const value_type& value)
    {
        return emplace_hint_key_args(hint, Traits::key(value), Traits::value(value));
    }
    iterator insert(iterator hint, value_type&& value)
    {
        return emplace_hint_key_args(hint, Traits::key(value), Traits::value(std::move(value)));
    }
    template <typename Arg> iterator insert(iterator hint, Arg&& value)
    {
        return emplace_hint_key_args(hint, Traits::key(value), Traits::value(std::move(value)));
    }

    template <typename InputIterator> void insert(InputIterator first, InputIterator last)
    {
        for (; first != last; ++first)
            insert(*first);
    }

    size_type count(fast_key_type key) const noexcept;

    size_type erase(fast_key_type key);
    iterator erase(iterator pos) { return pos != end() ? std::next(internal_erase(pos)) : pos; }

    // Function iterator erase(iterator first, iterator last) is not provided
    // because the ART iterators can invalidate, and there is no simple solution
    // for maintaining the last iterator un-invalidated in order to check the end
    // condition.

    void swap(self_t& other) noexcept;
    void clear();

    // Stats

    // Return current memory use by tree nodes in bytes.
    [[nodiscard]] constexpr size_type current_memory_use() const noexcept
    {
        return memory_use<leaf_type>() + memory_use<inode_4>() + memory_use<inode_16>() +
               memory_use<inode_64>() + memory_use<inode_256>();
    }

    [[nodiscard]] constexpr size_type leaf_count() const noexcept { return get_count<leaf_type>(); }
    [[nodiscard]] constexpr size_type inode4_count() const noexcept { return get_count<inode_4>(); }
    [[nodiscard]] constexpr size_type inode16_count() const noexcept
    {
        return get_count<inode_16>();
    }
    [[nodiscard]] constexpr size_type inode64_count() const noexcept
    {
        return get_count<inode_64>();
    }
    [[nodiscard]] constexpr size_type inode256_count() const noexcept
    {
        return get_count<inode_256>();
    }

    // Debugging
    void dump(std::ostream& os) const;

protected:
    template <typename... Args>
    [[nodiscard]] auto emplace_key_args(fast_key_type key, Args&&... args)
    {
        return emplace_key_args(multi_container(), bitwise_key(key), std::forward<Args>(args)...);
    }

    template <typename... Args>
    [[nodiscard]] iterator emplace_hint_key_args(iterator hint, fast_key_type key, Args&&... args);

private:
    using key_size_type = typename bitwise_key::size_type;
    using bitwise_key_prefix = std::pair<bitwise_key, key_size_type>;

    [[nodiscard]] const_iterator internal_locate(bitwise_key_prefix& key) const noexcept;
    [[nodiscard]] const_iterator internal_locate(bitwise_key_prefix&& key) const noexcept
    {
        return internal_locate(key);
    }

    [[nodiscard]] const_iterator internal_find(fast_key_type key) const noexcept;

    [[nodiscard]] const_iterator forward_step(const_iterator pos) const noexcept;
    template <typename Filter>
    [[nodiscard]] const_iterator internal_bound(bitwise_key key, Filter filter) const noexcept;
    [[nodiscard]] const_iterator internal_lower_bound(fast_key_type key) const noexcept;
    [[nodiscard]] const_iterator internal_upper_bound(fast_key_type key) const noexcept;

    template <typename... Args>
    [[nodiscard]] iterator internal_emplace(iterator hint, bitwise_key bitk,
                                            const bitwise_key_prefix& key, Args&&... args);

    template <typename... Args>
    [[nodiscard]] iterator emplace_key_args(std::true_type, bitwise_key key, Args&&... args);

    template <typename... Args>
    [[nodiscard]] std::pair<iterator, bool> emplace_key_args(std::false_type, bitwise_key key,
                                                             Args&&... args);

    iterator internal_erase(iterator pos);

private:
    using db_allocator_traits = std::allocator_traits<allocator_type>;

    [[nodiscard]] allocator_type& allocator() noexcept
    {
        return *static_cast<allocator_type*>(&tree);
    }

    template <typename Node> [[nodiscard]] constexpr size_type get_count() const noexcept
    {
        return std::get<counter<Node>>(count_).instances;
    }

    template <typename Node> [[nodiscard]] constexpr size_type memory_use() const noexcept
    {
        return sizeof(Node) * get_count<Node>();
    }

    template <typename Node> [[nodiscard]] constexpr size_type& count() noexcept
    {
        return std::get<counter<Node>>(count_).instances;
    }

    template <typename Node, typename... Args>
    unique_node_ptr<Node, self_t> make_node_ptr(Args&&... args);

    template <typename Node> void deallocate_node(Node* node) noexcept;
    void deallocate(inode_4* node) noexcept { deallocate_node(node); }
    void deallocate(inode_16* node) noexcept { deallocate_node(node); }
    void deallocate(inode_64* node) noexcept { deallocate_node(node); }
    void deallocate(inode_256* node) noexcept { deallocate_node(node); }
    void deallocate(leaf_type* leaf) noexcept(std::is_nothrow_destructible<mapped_type>::value);
    void deallocate(node_ptr node) noexcept(std::is_nothrow_destructible<mapped_type>::value);

    template <typename Node>
    void deallocate_subtree(Node* node) noexcept(std::is_nothrow_destructible<mapped_type>::value);

    void assign_to_parent(iterator hint, node_ptr child) noexcept;

    template <typename INode, typename NodePtr, typename SizeType>
    iterator create_inode(iterator hint, bitwise_key prefix, NodePtr pdst, leaf_unique_ptr leaf,
                          SizeType key);

    template <typename Source>
    iterator grow_node(iterator hint, leaf_unique_ptr leaf, std::uint8_t key_byte);

    // Functions to convert a node to a smaller type of node. This allows to fully
    // generalize the shrinking routine without stooping to strange hacks.
    template <typename Node>
    [[nodiscard]] unique_node_ptr<typename Node::smaller_inode_type, self_t> make_smaller_node(
        Node& src, iterator& pos)
    {
        auto inode = make_node_ptr<typename Node::smaller_inode_type>(src, pos.index());
        pos.parent_ = inode->tagged_self();
        pos.position = inode->index();
        return inode;
    }

    [[nodiscard]] static node_ptr make_smaller_node(inode_4& src, iterator& pos) noexcept
    {
        node_ptr last_child = src.leave_last_child(pos.index());
        if (src.parent() != nullptr) {
            pos.parent_ = src.parent();
            pos.position = pos.index() + src.index();
        } else {
            if (pos.index() == 0 && last_child.tag() != node_type::LEAF) {
                // Jump to the last child's tree
                pos.parent_ = last_child;
            } else {
                pos.node_ = last_child;
                pos.parent_ = nullptr;
            }
        }
        return last_child;
    }

    template <typename Node>
    [[nodiscard]] static node_ptr make_tagged_ptr(unique_node_ptr<Node, self_t> node) noexcept
    {
        return node_ptr(node.release(), Node::static_type());
    }

    [[nodiscard]] static constexpr node_ptr make_tagged_ptr(node_ptr p) noexcept { return p; }

    template <typename Source> iterator shrink_node(iterator pos);

private:
    // A helper struct to get the empty base class optimization for 0 size allocators.
    // In C++20 parlance that would be [[no_unique_address]] for the allocator.
    struct compress_empty_base : public allocator_type {
        compress_empty_base() = default;
        compress_empty_base(const allocator_type& alloc)
            : allocator_type(alloc)
        {
        }
        compress_empty_base(compress_empty_base&& rhs) = default;
        node_ptr root{};
    } tree;

    template <typename T> struct counter {
        size_type instances;
    };
    using node_stats_type = std::tuple<counter<leaf_type>, counter<inode_4>, counter<inode_16>,
                                       counter<inode_64>, counter<inode_256>>;
    node_stats_type count_{};
};

} // namespace detail
} // namespace art

#include "art_container_impl.h"

#endif // ART_DETAIL_ART_CONTAINER_HEADER_INCLUDED
