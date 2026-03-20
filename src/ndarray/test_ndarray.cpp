#include "ndarray.hpp"
#include <gtest/gtest.h>
// #include <vector>

TEST(NDArrayTest, ConstructWithShape)
{
    using Array3D = NDArray<int, 3, int>;
    Array3D::ShapeType shape = {2, 3, 4};
    Array3D arr(shape);

    EXPECT_EQ(arr.size(), 24);
    EXPECT_EQ(arr.shape()[0], 2);
    EXPECT_EQ(arr.shape()[1], 3);
    EXPECT_EQ(arr.shape()[2], 4);

    // 检查默认值
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 4; ++k)
            {
                // macro expand according comma, so parentheses are needed
                EXPECT_EQ((arr[i, j, k]), 0);
            }
}

TEST(NDArrayTest, ConstructWithPointer)
{
    using Array2D = NDArray<double, 2, int>;
    double data[6] = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
    Array2D::ShapeType shape = {2, 3};
    Array2D arr(shape, data);

    EXPECT_EQ(arr.size(), 6);
    EXPECT_EQ(arr.shape()[0], 2);
    EXPECT_EQ(arr.shape()[1], 3);

    // 检查数据是否正确引用
    auto x = arr[0, 0];
    EXPECT_EQ((arr[0, 0]), 1.1);
    x = arr[0, 1];
    EXPECT_EQ(x, 2.2);
    x = arr[1, 2];
    EXPECT_EQ(x, 6.6);
}

TEST(NDArrayTest, Slice)
{
    int data[24] = {};
    for (int i = 0; i < 24; ++i)
        data[i] = i;

    NDArray<int, 2, int> arr({4, 6}, data);

    auto slice = arr[Slice<int>{1, 3}, Slice<int>(2, 5)];
    EXPECT_EQ(slice.size(), 6);
    EXPECT_EQ(slice.shape()[0], 2);
    EXPECT_EQ(slice.shape()[1], 3);
    for (int i = 0; i < slice.shape()[0]; ++i)
    {
        for (int j = 0; j < slice.shape()[1]; ++j)
        {
            EXPECT_EQ((slice[i, j]), (arr[i + 1, j + 2]));
        }
    }

    auto slice2 = arr[1, Slice<int>(2, 6, 2)];
    EXPECT_EQ(slice2.size(), 2);
    EXPECT_EQ(slice2.shape()[0], 2);
    EXPECT_EQ((slice2[0]), (arr[1, 2]));
    EXPECT_EQ((slice2[1]), (arr[1, 4]));

    auto slice3 = arr[Slice<int>(3, 1, -1), Slice<int>(0, 4, 2)];
    EXPECT_EQ(slice3.size(), 4);
    EXPECT_EQ(slice3.shape()[0], 2);
    EXPECT_EQ(slice3.shape()[1], 2);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}