#pragma once
#ifndef DFT_ATOM_H
#define DFT_ATOM_H

#include <array>
#include <string>

namespace dft
{

class Atom
{
  public:
    using Vec3 = std::array<double, 3>;

    Atom() = default;

    Atom(const std::string &symbol, const Vec3 &position, double charge = 0.0)
        : symbol_(symbol), position_(position), charge_(charge)
    {
    }

    // 原子类型
    const std::string &symbol() const { return symbol_; }
    void set_symbol(const std::string &s) { symbol_ = s; }

    // 原子位置
    const Vec3 &position() const { return position_; }
    void set_position(const Vec3 &p) { position_ = p; }

    // 原子电荷（核电荷 / 有效电荷 / Bader 电荷由上层解释）
    double charge() const { return charge_; }
    void set_charge(double q) { charge_ = q; }

  private:
    std::string symbol_; // 元素符号
    Vec3 position_;      // (x,y,z)
    double charge_{0.0}; // 电荷
};

} // namespace dft

#endif // DFT_ATOM_H
