#ifndef CPPFLOW_STUB_HPP
#define CPPFLOW_STUB_HPP

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace cppflow
{

class tensor
{
  public:
    tensor(int32_t) {}
    tensor(std::vector<int32_t>, std::vector<int32_t>) {}
    tensor(std::vector<double>, std::vector<int32_t>) {}
};

class model
{
  public:
    model(std::string) {};
    std::vector<tensor> operator()(std::vector<std::tuple<std::string, tensor>>,
                                   std::vector<std::string>)
    {
        return std::vector<tensor>();
    }
};

} // namespace cppflow

#endif
