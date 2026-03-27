#include "ndarray/ndarray.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
bool is_equal(double lhs, double rhs)
{
    return std::abs(lhs - rhs) < 1e-12;
}
} // namespace

int main()
{
    {
        using Array3D = NDArray<int, 3, int>;
        Array3D::ShapeType shape = {2, 3, 4};
        Array3D array(shape);

        if (array.size() != 24)
        {
            std::cerr << "Unexpected 3D array size: " << array.size() << '\n';
            return 1;
        }

        if (array.shape()[0] != 2 || array.shape()[1] != 3 || array.shape()[2] != 4)
        {
            std::cerr << "Unexpected 3D array shape" << '\n';
            return 1;
        }

        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 4; ++k)
                {
                    if (array[i, j, k] != 0)
                    {
                        std::cerr << "Default-initialized element is not zero" << '\n';
                        return 1;
                    }
                }
            }
        }
    }

    {
        auto default_heap_array = NDArray<int, 2, int>::with_static_shape<2, 3>();
        if (default_heap_array.uses_stack_storage())
        {
            std::cerr << "Default max_stack_elements=64 should force heap storage" << '\n';
            return 1;
        }

        auto explicit_heap_array = NDArray<int, 2, int, 64>::with_static_shape<2, 3>();
        if (explicit_heap_array.uses_stack_storage())
        {
            std::cerr << "Explicit max_stack_elements=64 should force heap storage" << '\n';
            return 1;
        }

        auto custom_stack_array = NDArray<int, 2, int, 8>::with_static_shape<2, 4>();
        if (!custom_stack_array.uses_stack_storage())
        {
            std::cerr << "Custom stack-buffer limit was not applied" << '\n';
            return 1;
        }

        bool saw_stack_capacity_error = false;
        try
        {
            NDArray<int, 2, int, 8> oversized_stack_array({3, 3});
        }
        catch (const std::length_error &)
        {
            saw_stack_capacity_error = true;
        }

        if (!saw_stack_capacity_error)
        {
            std::cerr << "Stack-selected array did not reject oversized runtime shape" << '\n';
            return 1;
        }

        NDArray<int, 2, int, 0> zero_stack_limit_array({1, 1});
        if (zero_stack_limit_array.uses_stack_storage())
        {
            std::cerr << "max_stack_elements=0 should force heap storage" << '\n';
            return 1;
        }
    }

    {
        NDArray<int, 2, int, 8> stack_source({2, 4});
        stack_source[0, 0] = 7;
        NDArray<int, 2, int, 8> stack_copy(stack_source);
        if (!stack_copy.uses_stack_storage() || stack_copy[0, 0] != 7)
        {
            std::cerr << "Stack-backed copy did not preserve stack storage" << '\n';
            return 1;
        }

        NDArray<int, 2, int, 8> stack_moved(std::move(stack_source));
        if (!stack_moved.uses_stack_storage() || stack_moved[0, 0] != 7)
        {
            std::cerr << "Stack-backed move did not preserve stack storage" << '\n';
            return 1;
        }

        NDArray<int, 2, int> heap_source({2, 3});
        heap_source[1, 2] = 19;
        NDArray<int, 2, int> heap_copy(heap_source);
        if (heap_copy.uses_stack_storage() || heap_copy[1, 2] != 19)
        {
            std::cerr << "Heap-backed copy did not preserve heap storage" << '\n';
            return 1;
        }

        NDArray<int, 2, int> heap_moved(std::move(heap_source));
        if (heap_moved.uses_stack_storage() || heap_moved[1, 2] != 19)
        {
            std::cerr << "Heap-backed move did not preserve heap storage" << '\n';
            return 1;
        }
    }

    {
        using Array2D = NDArray<double, 2, int>;
        double data[6] = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
        Array2D::ShapeType shape = {2, 3};
        Array2D array(shape, data);

        if (array.size() != 6 || array.shape()[0] != 2 || array.shape()[1] != 3)
        {
            std::cerr << "Unexpected pointer-backed array shape" << '\n';
            return 1;
        }

        if (!is_equal(array[0, 0], 1.1) || !is_equal(array[0, 1], 2.2) || !is_equal(array[1, 2], 6.6))
        {
            std::cerr << "Pointer-backed array did not expose expected values" << '\n';
            return 1;
        }
    }

    {
        bool saw_negative_shape = false;
        try
        {
            NDArray<int, 2, int> invalid_shape_array({-1, 3});
        }
        catch (const std::invalid_argument &)
        {
            saw_negative_shape = true;
        }

        if (!saw_negative_shape)
        {
            std::cerr << "Negative shape was not rejected" << '\n';
            return 1;
        }
    }

    {
        bool saw_null_pointer = false;
        try
        {
            NDArray<int, 2, int> invalid_view({2, 2}, nullptr);
        }
        catch (const std::invalid_argument &)
        {
            saw_null_pointer = true;
        }

        if (!saw_null_pointer)
        {
            std::cerr << "Null data pointer for non-empty shape was not rejected" << '\n';
            return 1;
        }
    }

#ifdef TEST_COMPILE
    {
        bool saw_shape_overflow = false;
        try
        {
            NDArray<int, 2, int> overflow_shape_array({46341, 46341});
        }
        catch (const std::overflow_error &)
        {
            saw_shape_overflow = true;
        }

        if (!saw_shape_overflow)
        {
            std::cerr << "Overflowing shape was not rejected in test mode" << '\n';
            return 1;
        }
    }
#endif

    {
        int data[24] = {};
        for (int i = 0; i < 24; ++i)
        {
            data[i] = i;
        }

        NDArray<int, 2, int> array({4, 6}, data);
        array[2, 3] = -17;
        if (array[0, 0] != 0 || array[1, 4] != 10 || array[2, 3] != -17)
        {
            std::cerr << "Basic multidimensional indexing failed" << '\n';
            return 1;
        }

        const NDArray<int, 2, int> &const_array = array;
        if (const_array[2, 3] != -17)
        {
            std::cerr << "Const multidimensional indexing failed" << '\n';
            return 1;
        }

        bool saw_out_of_bounds = false;
        try
        {
            (void)array[4, 0];
        }
        catch (const std::out_of_range &)
        {
            saw_out_of_bounds = true;
        }

        if (!saw_out_of_bounds)
        {
            std::cerr << "Out-of-bounds indexing was not detected in test mode" << '\n';
            return 1;
        }
    }

    return 0;
}
