#include <iostream>

template <typename T, T t>
T  tfunc()
{
    return t + 10;
}

template <typename T, T t>
constexpr T  func()
{
    return tfunc<T, t>();
}

#define FUNC(a)  func<decltype(a),a>()

int main()
{
    std::cout << FUNC(10) << std::endl;
}