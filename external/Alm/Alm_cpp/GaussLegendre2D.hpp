#pragma once

/** @file GaussLegendre.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 08/04/17
  */

#include<cstddef>
#include<array>
#include "GaussLegendre1D.hpp"

namespace _2D {
namespace GQ {

template<typename T, std::size_t Order>
class GaussLegendreQuadrature
{
  public:
    _1D::GQ::GaussLegendreQuadrature<T,Order> _1dInt;

    GaussLegendreQuadrature() = default;

    // This version will integrate a callable between four points + 4 additional parameters for the function
    template<typename F, typename X, typename Y, typename Z , typename ZZ, typename ZZZ, typename ZZZZ>
    T operator()( F f, X a, X b, Y c, Y d , Z p0, ZZ p1, ZZZ p2, ZZZZ p3) const; 
  protected:
};


template<typename T, std::size_t Order>
template<typename F, typename X, typename Y, typename Z, typename ZZ, typename ZZZ, typename ZZZZ>
T GaussLegendreQuadrature<T,Order>::operator()(F f, X a, X b, Y c, Y d, Z p0, ZZ p1, ZZZ p2, ZZZZ p3) const
{
  // A 2D integral I = \int \int f(x,y) dx dy
  // can be written as two 1D integrals
  //
  // g(x) = \int f(x,y) dy
  // I = \int g(x) dx
  //
  // first, integrate over y at each of the required points (i.e. create g(x_i))

  X apb = (b + a)/2;
  X amb = (b - a)/2;
  std::array<T, Order> sums;

  //std::cout << "v =" << v << std::endl;
  //std::cout << "w =" << w << std::endl;
  #pragma parallel for
  for(std::size_t i = 0; i < Order; i++)
    sums[i] = _1dInt( [&](Y y){ return f(apb + amb*_1dInt.getX()[i], y, p0 ,p1, p2, p3); }, c, d);

  // now integrate over x
  T sum = 0;
  for(std::size_t i = 0; i < Order; i++)
    sum += _1dInt.getW()[i]*sums[i];
  sum *= amb;

  return sum;
}

}

}
