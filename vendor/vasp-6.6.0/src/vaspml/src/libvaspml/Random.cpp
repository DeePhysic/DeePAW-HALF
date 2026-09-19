#include "Random.hpp"

#include <algorithm>  // for shuffle
#include <iostream>   // for basic_ostream, operator<<, basic_ostream::operator<<, cerr, cout
#include <numeric>    // for iota

using namespace vaspml;

UniqueRandomizer::UniqueRandomizer() : generator(std::random_device{}())
{
    minVal = 0;
    maxVal = 0;
}

UniqueRandomizer::UniqueRandomizer(const Size seed) : generator(seed)
{
    minVal = 0;
    maxVal = 0;
}

UniqueRandomizer::UniqueRandomizer(const Vec1Int& seed)
{
    minVal = 0;
    maxVal = 0;
    std::seed_seq seq(seed.cbegin(), seed.cend());
    generator.seed(seq);
}

void UniqueRandomizer::init(const Size seed)
{
    minVal = 0;
    maxVal = 0;
    generator.seed(seed);

    return;
}

void UniqueRandomizer::init(const Vec1Int& seed)
{
    minVal = 0;
    maxVal = 0;
    std::seed_seq seq(seed.cbegin(), seed.cend());
    generator.seed(seq);

    return;
}

UniqueRandomizer::UniqueRandomizer(const Int min, const Int max) : generator(std::random_device{}())
{
    initializeNumbers(min, max);
}

UniqueRandomizer::UniqueRandomizer(const Int min, const Int max, const Size seed) : generator(seed)
{
    initializeNumbers(min, max);
}

void UniqueRandomizer::reinitializeRange(const Int newMin, const Int newMax)
{
    if (newMin > newMax)
    {
        std::cerr << "Warning: new_min is greater than new_max. Range not updated.\n";
        return;
    }
    std::cout << "Range successfully reinitialized from [" << minVal << ", " << maxVal << ") to ["
              << newMin << ", " << newMax << ").\n";
    initializeNumbers(newMin, newMax);

    return;
}

Vec1Int UniqueRandomizer::draw(const Int count)
{
    if ((Size)count >= numbers.size())
    {
        std::cerr << "Error: Draw count (" << count << ") exceeds current range size ("
                  << numbers.size() << ").\n";
        return numbers;
    }
    // 1. Shuffle the entire vector to randomize the order
    std::shuffle(numbers.begin(), numbers.end(), generator);
    // 2. Select the required count of numbers
    Vec1Int uniqueNumbers(numbers.cbegin(), numbers.cbegin() + count);

    return uniqueNumbers;
}

void UniqueRandomizer::initializeNumbers(const Int min, const Int max)
{
    minVal = min;
    maxVal = max;
    int range_size = maxVal - minVal;

    // Clear and resize the vector
    numbers.clear();
    numbers.resize(range_size);
    // Use std::iota to efficiently fill the vector with the sequence [min_val, ..., max_val]
    std::iota(numbers.begin(), numbers.end(), minVal);

    return;
}

UniformRandomIntGenerator::UniformRandomIntGenerator()
{
    std::random_device rd;
    generator.seed(rd());
}

UniformRandomIntGenerator::UniformRandomIntGenerator(const Size seed)
{
    generator.seed(seed);
}

UniformRandomIntGenerator::UniformRandomIntGenerator(const Vec1Int& seed)
{
    std::seed_seq seq(seed.cbegin(), seed.cend());
    generator.seed(seq);
}

void UniformRandomIntGenerator::init(const Size seed)
{
    generator.seed(seed);

    return;
}

void UniformRandomIntGenerator::init(const Vec1Int& seed)
{
    std::seed_seq seq(seed.cbegin(), seed.cend());
    generator.seed(seq);

    return;
}

Int UniformRandomIntGenerator::uniformInt(Int low, Int high)
{
    std::uniform_int_distribution<> distrib(low, high);

    return distrib(generator);
}
