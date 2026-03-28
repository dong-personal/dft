#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "mdspan/mdspan/mdspan.hpp"

#ifndef RankType
#define RankType int
#endif // RankType

// Lightweight local array wrapper for small in-process data.
//
// Notice:
// - This class only supports basic multidimensional indexing.
// - Slice, subview, and MPI/distributed-data use cases are intentionally not supported.
// - Use MFEM containers for data that may need MPI communication or distributed ownership.
// - This wrapper assumes contiguous row-major storage for its owning-storage path.
// - The stack-buffer capacity is a compile-time template parameter.
template <typename T, RankType rank, typename Index, std::size_t max_stack_elements = 0>
class NDArray
{
  public:
    using ShapeType = std::array<Index, rank>;

    template <typename, RankType, typename, std::size_t>
    friend class NDArray;

  private:
    using MdspanType = Kokkos::mdspan<T, Kokkos::dextents<Index, rank>, Kokkos::layout_stride>;

    static constexpr std::size_t kMaxStackElements = max_stack_elements;
    static constexpr bool kUsesStackStorage = kMaxStackElements > 0 && kMaxStackElements < 64;

    struct HeapStorage
    {
        std::vector<T> heap_data;

        T *data_pointer() { return heap_data.data(); }

        const T *data_pointer() const { return heap_data.data(); }

        void assign(std::size_t element_count) { heap_data.assign(element_count, T{}); }
    };

    struct StackStorage
    {
        std::array<T, kMaxStackElements> stack_data{};

        T *data_pointer() { return stack_data.data(); }

        const T *data_pointer() const { return stack_data.data(); }

        void assign(std::size_t element_count) { std::fill_n(stack_data.begin(), element_count, T{}); }
    };

    using OwningStorage = std::conditional_t<kUsesStackStorage, StackStorage, HeapStorage>;

    std::unique_ptr<OwningStorage> m_storage;
    MdspanType m_mdspan;

    static constexpr Index size_from_shape(const ShapeType &shape)
    {
        Index total = 1;
        for (const Index dimension : shape)
        {
            total *= dimension;
        }
        return total;
    }

    static constexpr ShapeType stride_from_shape(const ShapeType &shape)
    {
        ShapeType stride{};
        for (RankType i = 0; i < rank; ++i)
        {
            stride[i] = 1;
            for (RankType j = i + 1; j < rank; ++j)
            {
                stride[i] *= shape[j];
            }
        }
        return stride;
    }

    static constexpr bool has_negative_extent(const ShapeType &shape)
    {
        if constexpr (std::is_signed_v<Index>)
        {
            for (const Index dimension : shape)
            {
                if (dimension < 0)
                {
                    return true;
                }
            }
        }
        return false;
    }

#ifdef TEST_COMPILE
    static constexpr bool size_overflows_index(const ShapeType &shape)
    {
        std::size_t total = 1;
        for (const Index dimension : shape)
        {
            const std::size_t extent = static_cast<std::size_t>(dimension);
            if (extent != 0 && total > std::numeric_limits<std::size_t>::max() / extent)
            {
                return true;
            }
            total *= extent;
        }

        return total > static_cast<std::size_t>(std::numeric_limits<Index>::max());
    }
#endif

    static void validate_shape(const ShapeType &shape)
    {
#ifdef TEST_COMPILE
        if (has_negative_extent(shape))
        {
            throw std::invalid_argument("NDArray shape dimensions must be non-negative");
        }

        if (size_overflows_index(shape))
        {
            throw std::overflow_error("NDArray shape size overflows index type");
        }
#endif
    }

    static std::size_t element_count_from_shape(const ShapeType &shape)
    {
        validate_shape(shape);

        std::size_t total = 1;
        for (const Index dimension : shape)
        {
            total *= static_cast<std::size_t>(dimension);
        }
        return total;
    }

    static MdspanType make_mdspan(T *data, const ShapeType &shape)
    {
        return MdspanType(
            data, Kokkos::layout_stride::mapping<Kokkos::dextents<Index, rank>>(shape, stride_from_shape(shape)));
    }

    template <typename... Indices>
    static constexpr bool kHasValidIndices =
        sizeof...(Indices) == rank && (std::is_integral_v<std::decay_t<Indices>> && ...);

    template <typename... Indices>
    bool indices_in_bounds(Indices... indices) const
        requires(kHasValidIndices<Indices...>)
    {
        const std::array<Index, rank> index_array = {static_cast<Index>(indices)...};
        for (RankType i = 0; i < rank; ++i)
        {
            if (index_array[i] < 0 || index_array[i] >= m_mdspan.extent(i))
            {
                return false;
            }
        }
        return true;
    }

#ifdef TEST_COMPILE
    template <typename... Indices>
    void check_indices_in_bounds(Indices... indices) const
        requires(kHasValidIndices<Indices...>)
    {
        if (!indices_in_bounds(indices...))
        {
            throw std::out_of_range("NDArray index out of bounds");
        }
    }
#endif

    T *owning_data_pointer() { return m_storage ? m_storage->data_pointer() : nullptr; }

    const T *owning_data_pointer() const { return m_storage ? m_storage->data_pointer() : nullptr; }

    void bind_mdspan(T *data, const ShapeType &shape)
    {
        const std::size_t element_count = element_count_from_shape(shape);
#ifdef TEST_COMPILE
        if (element_count > 0 && data == nullptr)
        {
            throw std::invalid_argument("NDArray data pointer must not be null for a non-empty shape");
        }
#endif

        m_mdspan = make_mdspan(data, shape);
    }

    void assign_owned_storage(const ShapeType &shape)
    {
        const std::size_t element_count = element_count_from_shape(shape);
#ifdef TEST_COMPILE
        if constexpr (kUsesStackStorage)
        {
            if (element_count > kMaxStackElements)
            {
                throw std::length_error("NDArray shape exceeds compile-time stack capacity");
            }
        }
#endif

        m_storage = std::make_unique<OwningStorage>();
        m_storage->assign(element_count);
        bind_mdspan(owning_data_pointer(), shape);
    }

    void copy_from_other(const NDArray &other)
    {
        assign_owned_storage(other.shape());
        std::copy_n(other.data_pointer(), other.size(), owning_data_pointer());
    }

    void move_from_other(NDArray &&other)
    {
        if (!other.data_pointer())
        {
            m_mdspan = MdspanType();
            m_storage.reset();
            return;
        }

        if (!other.m_storage)
        {
            m_mdspan = other.m_mdspan;
            m_storage.reset();
            return;
        }

        const ShapeType other_shape = other.shape();
        if constexpr (kUsesStackStorage)
        {
            m_storage = std::make_unique<OwningStorage>();
            m_storage->assign(other.size());
            std::copy_n(other.owning_data_pointer(), other.size(), m_storage->data_pointer());
            bind_mdspan(m_storage->data_pointer(), other_shape);
            return;
        }

        m_storage = std::move(other.m_storage);
        bind_mdspan(m_storage->data_pointer(), other_shape);
    }

  public:
    NDArray() = default;

    template <Index... static_shape>
        requires(sizeof...(static_shape) == rank)
    static auto with_static_shape()
    {
        constexpr ShapeType shape = {static_shape...};
        static_assert(!has_negative_extent(shape), "Negative static shape is not allowed");
#ifdef TEST_COMPILE
        if constexpr (!has_negative_extent(shape))
        {
            static_assert(!size_overflows_index(shape), "Static shape size overflows index type");
        }
#endif
        constexpr Index static_size = size_from_shape(shape);
        if constexpr (kUsesStackStorage)
        {
            static_assert(static_size <= kMaxStackElements, "Static shape exceeds stack storage capacity");
        }
        NDArray<T, rank, Index, kMaxStackElements> array;
        array.assign_owned_storage(shape);
        return array;
    }

    // Owning constructor for local contiguous storage.
    explicit NDArray(const ShapeType &shape) { assign_owned_storage(shape); }

    // Non-owning view constructor over external contiguous storage.
    NDArray(const ShapeType &shape, T *data) { bind_mdspan(data, shape); }

    explicit NDArray(MdspanType mdspan) : m_mdspan(mdspan) {}

    NDArray(const NDArray &other) { copy_from_other(other); }

    NDArray(NDArray &&other) noexcept { move_from_other(std::move(other)); }

    NDArray &operator=(const NDArray &other)
    {
        if (this != &other)
        {
            copy_from_other(other);
        }
        return *this;
    }

    NDArray &operator=(NDArray &&other) noexcept
    {
        if (this != &other)
        {
            move_from_other(std::move(other));
        }
        return *this;
    }

    Index size() const { return m_mdspan.size(); }

    ShapeType shape() const
    {
        ShapeType shape{};
        for (RankType i = 0; i < rank; ++i)
        {
            shape[i] = m_mdspan.extent(i);
        }
        return shape;
    }

    T *data_pointer() { return m_mdspan.data_handle(); }

    const T *data_pointer() const { return m_mdspan.data_handle(); }

    bool uses_stack_storage() const { return m_storage && kUsesStackStorage; }

    template <typename... Indices>
        requires(kHasValidIndices<Indices...>)
    decltype(auto) operator[](Indices... indices)
    {
#ifdef TEST_COMPILE
        check_indices_in_bounds(indices...);
#endif
        // Only full-rank integer indexing is allowed.
        return m_mdspan[indices...];
    }

    template <typename... Indices>
        requires(kHasValidIndices<Indices...>)
    decltype(auto) operator[](Indices... indices) const
    {
#ifdef TEST_COMPILE
        check_indices_in_bounds(indices...);
#endif
        // Only full-rank integer indexing is allowed.
        return m_mdspan[indices...];
    }
};

// Hot fixed-size instantiations used by Atom/Structure are emitted once in a
// source file to reduce repeated template code generation across translation
// units.
#ifndef TEST_COMPILE
extern template class NDArray<double, 1, int, 3>;
extern template class NDArray<double, 2, int, 9>;
#endif

#endif // NDARRAY_HPP
