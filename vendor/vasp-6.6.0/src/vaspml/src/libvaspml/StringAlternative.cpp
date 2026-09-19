#include "StringAlternative.hpp"

#include "Record.hpp"
#include "Tutor.hpp"

#include <algorithm> // for find
#include <iterator>  // for distance, __iterator_category

using namespace vaspml;

MapString vaspml::dataMapStringAlternative()
{
    return MapString{
        {"numUnique", "Int"       },
        {"astring",   "Vec1String"},
        {"index",     "Vec1Int"   }
    };
}

StringAlternative::StringAlternative(ShRec record) :
    data(assignOrMakeRecord(record, dataMapStringAlternative())),
    numUnique(data->get<Int>("numUnique")),
    astring(data->get<Vec1String>("astring")),
    index(data->get<Vec1Int>("index"))
{}

String StringAlternative::operator()(const String& arg) const
{
    // Search for string in list of all alternatives.
    Size i = std::distance(astring.begin(), std::find(astring.begin(), astring.end(), arg));
    if (i >= astring.size()) return arg;
    // If already main alternative, return with input string.
    else if (index[i] < 0) return arg;
    // Else, search for negative index in list of all indices.
    else return searchMain(i, arg);
}

String StringAlternative::operator[](const String& arg) const
{
    // Search for string in list of all alternatives.
    Size i = std::distance(astring.begin(), std::find(astring.begin(), astring.end(), arg));
    if (i >= astring.size())
    {
        global::tutor.bug("Cannot find string alternative to \"" + arg + "\".");
    }
    // If already main alternative, return with input string.
    if (index[i] < 0) return arg;
    // Else, search for negative index in list of all indices.
    else return searchMain(i, arg);
}

void StringAlternative::add(const Vec2String& alternatives)
{
    for (const auto& a : alternatives) add(a);

    return;
}

void StringAlternative::add(const Vec1String& alternatives)
{
    if (alternatives.empty()) return;

    numUnique++;
    for (Int i = 0; i < static_cast<Int>(alternatives.size()); ++i)
    {
        if (alternatives[i].empty()) continue;
        astring.push_back(alternatives[i]);
        if (i == 0) index.push_back(-numUnique);
        else index.push_back(numUnique);
    }

    return;
}

String StringAlternative::searchMain(Size i, String arg) const
{
    // Else, search for negative index in list of all indices.
    Int  n = -index[i];
    Size j = std::distance(index.begin(), std::find(index.begin(), index.end(), n));
    if (j >= astring.size())
    {
        global::tutor.bug("Invalid index encountered for string alternative to \"" + arg + "\".");
    }

    return astring[j];
}
