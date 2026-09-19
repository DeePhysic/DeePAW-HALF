#include "Tutor.hpp"

#include <iostream>
#include <stdexcept>

using namespace vaspml;

void Tutor::warning(const String message) const
{
    std::cout << headerWarning;
    std::cout << message;
    std::cout << footerWarning;

    return;
}

void Tutor::error(const String message) const
{
    std::cout << headerError;
    std::cout << message;
    std::cout << footerError;
    throw std::runtime_error("");
}

void Tutor::bug(const String message) const
{
    std::cout << headerBug;
    std::cout << message;
    std::cout << footerBug;
    throw std::runtime_error("");
}

namespace global
{
Tutor tutor;
}
