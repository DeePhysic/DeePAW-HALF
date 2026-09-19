#ifndef BOOLFLAVOR_HPP
#define BOOLFLAVOR_HPP

#include "types.hpp"

namespace vaspml::text
{

class BoolFlavor
{
  public:
    BoolFlavor(Vec1String trues, Vec1String falses, bool fortranStyle);
    virtual bool isBool(const String& input) const;
    virtual bool convert(const String& input) const;

  private:
    Vec1String trues;
    Vec1String falses;
    bool       fortranStyle;
};

class BoolFlavorJson : public BoolFlavor
{
  public:
    BoolFlavorJson();
};

class BoolFlavorExtXyz : public BoolFlavor
{
  public:
    BoolFlavorExtXyz();
};

class BoolFlavorIncar : public BoolFlavor
{
  public:
    BoolFlavorIncar();
};

} //namespace vaspml::text

#endif
