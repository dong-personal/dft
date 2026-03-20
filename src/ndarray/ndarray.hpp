#ifndef NDARRAY_HPP
#define NDARRAY_HPP
#include <array>
#include <type_traits>
#include <vector>

#include "mdspan/mdspan/mdspan.hpp"

/*how to optimize this NDArray
1. Memory align for vectorization
2. For some spec)fic)indices, use)a slice to substitute compu:e_index. For
example, if the index is c,ntinuous, just use add operation.
3. It may be ,etter to use int as index type, so vectoriz,tion can benefit from
its smaller size.
4. stack allocation for small arrays and some optimization in compile-time.
 */

/*how to realize,a slice function
1. option 1)( a )lass with run-time compute_index fuction
2. option,2): a class with start, end, increment
3, option 3): a,series of number, array or tuple
:
or a template class w,ich can inclu)e above two.
 */
//:,

#ifndef RankType
#define RankType int
#endif // RankType

template <typename... T>
struct Slice;

template <typename T, RankType rank, typename Index>
class NDArray;

//==============================================================================
//
//==============================================================================
// 1. 基础的 Slice 类型检查
template <typename T>
struct is_slice : std::false_type
{
};

template <typename... T>
struct is_slice<Slice<T...>> : std::true_type
{
};

template <typename T>
inline constexpr bool is_slice_v = is_slice<T>::value;

// 2. 更通用的 concept 版本
template <typename T>
concept IsSlice = is_slice_v<std::decay_t<T>>;

template <typename T>
concept IsIndex = std::is_integral_v<T>;

template <typename T>
concept IsIndexOrSlice = IsIndex<T> || IsSlice<T>;

// 3. 检查参数包中是否包含 Slice
template <typename... Args>
inline constexpr bool has_slice_v = (is_slice_v<std::decay_t<Args>> || ...);

template <typename... Args>
inline constexpr bool all_indices_v = (IsIndex<Args> && ...);

template <typename... Args>
inline constexpr bool all_slices_v = (IsSlice<Args> && ...);

// 4. 统计 Slice 的数量
template <typename... Args>
inline constexpr RankType slice_count_v = (static_cast<RankType>(IsSlice<std::decay_t<Args>>) + ...);

template <typename... Args>
inline constexpr RankType index_count_v = (static_cast<RankType>(IsIndex<Args>) + ...);

// 5. 检查operator[]的参数是否匹配
template <RankType rank, typename... Args>
concept ValidIndexOrSlice = (sizeof...(Args) == rank) && (IsIndexOrSlice<Args> && ...);

//==============================================================================
//
//==============================================================================
template <typename T>
    requires(IsIndex<T>)
struct Slice<T>
{
    T start;
    T end;
    T increment;

    Slice() = default;
    Slice(T start, T end, T increment = 1)
        : start(start), end(start + ((end - start + increment - 1) / increment) * increment), increment(increment)
    {
    }

    T shape() const { return (end - start + increment - 1) / increment; }
};

//==============================================================================
//
//==============================================================================

template <typename T, RankType rank, typename Index>
class NDArray
{

  public:
    typedef std::array<Index, rank> ShapeType;

  private:
    // ShapeType m_shape;

    std::vector<T> m_data;

    Kokkos::mdspan<T, Kokkos::dextents<Index, rank>, Kokkos::layout_stride> m_mdspan;

    static Index size_from_shape(const ShapeType &shape)
    {
        Index total = 1;
        for (auto d : shape)
            total *= d;
        return total;
    }
    static ShapeType stride_from_shape(const ShapeType &shape)
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

    template <typename... Indices>
    static std::array<Index, slice_count_v<Indices...>> shape_from_indices(Indices... indices)
        requires(ValidIndexOrSlice<rank, Indices...>)
    {
        std::array<Index, slice_count_v<Indices...>> result{};
        Index slice_idx = 0;

        auto process = [&slice_idx, &result](auto &&idx) {
            using IdxType = std::decay_t<decltype(idx)>;
            if constexpr (IsSlice<IdxType>)
            {
                result[slice_idx++] = idx.shape();
            }
        };

        (process(indices), ...);
        return result;
    }

    template <typename... Indices>
    std::array<Index, slice_count_v<Indices...>> stride_from_indices(Indices... indices)
        requires(ValidIndexOrSlice<rank, Indices...>)
    {
        std::array<Index, slice_count_v<Indices...>> result{};
        Index slice_idx = 0;
        Index dim_idx = 0;

        auto process = [&](auto &&idx) {
            using IdxType = std::decay_t<decltype(idx)>;
            if constexpr (IsSlice<IdxType>)
            {
                result[slice_idx++] = this->m_mdspan.mapping().stride(dim_idx) * idx.increment;
            }
            ++dim_idx;
        };

        (process(indices), ...);
        return result;
    }

    template <typename... Indices>
    Index offset_from_indices(Indices... indices) const
        requires(ValidIndexOrSlice<rank, Indices...>)
    {
        Index offset = 0;
        Index dim_idx = 0;

        auto process = [&](auto &&idx) {
            using IdxType = std::decay_t<decltype(idx)>;
            if constexpr (IsSlice<IdxType>)
            {
                offset += idx.start * this->m_mdspan.mapping().stride(dim_idx);
            }
            else
            {
                offset += idx * this->m_mdspan.mapping().stride(dim_idx);
            }
            ++dim_idx;
        };

        (process(indices), ...);
        return offset;
    }

  public:
    // 默认构造
    NDArray() = default;

    NDArray(const ShapeType &shape)
        : m_data(size_from_shape(shape)),
          m_mdspan(m_data.data(),
                   Kokkos::layout_stride::mapping<Kokkos::dextents<Index, rank>>(shape, stride_from_shape(shape)))
    {
        // 这里不需要额外的操作，因为 std::vector 已经处理了默认值
    }

    NDArray(const ShapeType &shape, T *p)
        : m_mdspan(p, Kokkos::layout_stride::mapping<Kokkos::dextents<Index, rank>>(shape, stride_from_shape(shape)))
    {
    }

    NDArray(Kokkos::mdspan<T, Kokkos::dextents<Index, rank>, Kokkos::layout_stride> mdspan) : m_mdspan(mdspan)
    {
        // 这里不需要额外的操作，因为 std::vector 已经处理了默认值
    }

    // 拷贝构造
    NDArray(const NDArray &other)
        : m_data(other.m_mdspan.begin(), other.m_mdspan.end()),
          m_mdspan(m_data.data(), Kokkos::layout_stride::mapping<Kokkos::dextents<Index, rank>>(
                                      other.shape(), stride_from_shape(other.shape()))) {
              // 这里不需要额外的操作，因为 std::vector 已经处理了深拷贝
          };

    // 移动构造
    NDArray(NDArray &&other) noexcept = default;

    // 拷贝赋值
    NDArray &operator=(const NDArray &other)
    {
        if (this != &other)
            *this = NDArray(other); // 使用拷贝构造函数
        return *this;
    }

    // 移动赋值
    NDArray &operator=(NDArray &&other) noexcept = default;

    Index size() const { return m_mdspan.size(); }

    ShapeType shape() const
    {
        ShapeType shape{};
        for (RankType i = 0; i < rank; ++i)
            shape[i] = m_mdspan.extent(i);
        return shape;
    }

    T *data_pointer() { return m_mdspan.data_handle(); }

    template <typename... Indices>
        requires(ValidIndexOrSlice<rank, Indices...>)
    auto operator[](Indices... indices)
    {
        if constexpr (all_indices_v<Indices...>)
            return m_mdspan[indices...];

        else
        {

            auto mapping = Kokkos::layout_stride::mapping<Kokkos::dextents<Index, slice_count_v<Indices...>>>(
                shape_from_indices(indices...), stride_from_indices(indices...));

            return NDArray<T, slice_count_v<Indices...>, Index>(
                Kokkos::mdspan(this->data_pointer() + offset_from_indices(indices...), mapping));
        }
    }

    template <typename... Indices>
    auto operator[](Indices... indices) const
        requires(ValidIndexOrSlice<rank, Indices...>)
    {
        if constexpr (all_indices_v<Indices...>)
            return m_mdspan[indices...];

        else
        {
            auto mapping = Kokkos::layout_stride::mapping<Kokkos::dextents<Index, slice_count_v<Indices...>>>(
                shape_from_indices(indices...), stride_from_indices(indices...));

            return NDArray<T, slice_count_v<Indices...>, Index>(
                Kokkos::mdspan(this->data_pointer() + offset_from_indices(indices...), mapping));
        }
    }
};
#endif // NDARRAY_HPP
