#include "Timer.hpp"

#include "Record.hpp"
#include "Tutor.hpp"
#include "utils.hpp"

#include <algorithm>   // for find, max, min
#include <iostream>    // for basic_ostream, cout, operator<<
#include <iterator>    // for distance
#include <limits>      // for numeric_limits

using namespace vaspml;
using hrc = std::chrono::high_resolution_clock;

MapString vaspml::dataMapTimer()
{
    return MapString{
        {"ids",     "Vec1String"},
        {"tSum",    "Vec1Real"  },
        {"tMin",    "Vec1Real"  },
        {"tMax",    "Vec1Real"  },
        {"nStarts", "Vec1Real"  },
        {"status",  "Vec1Int"   }
    };
}

Timer::Timer(ShRec record) :
    data(assignOrMakeRecord(record, dataMapTimer())),
    ids(data->get<Vec1String>("ids")),
    tSum(data->get<Vec1Real>("tSum")),
    tMin(data->get<Vec1Real>("tMin")),
    tMax(data->get<Vec1Real>("tMax")),
    nStarts(data->get<Vec1Real>("nStarts")),
    status(data->get<Vec1Int>("status"))
{}

void Timer::start(const String& id)
{
    // Take time immediately, check if timer exists later.
    auto current = hrc::now();

    auto pos = std::find(ids.begin(), ids.end(), id);

    if (pos == ids.end())
    {
        ids.push_back(id);
        tSum.push_back(0.0);
        tMin.push_back(std::numeric_limits<Real>::max());
        tMax.push_back(std::numeric_limits<Real>::lowest());
        nStarts.push_back(1.0);
        status.push_back(1);
        tStart.push_back(current);
        tStop.push_back(current);
    }
    else
    {
        const Size i = std::distance(ids.begin(), pos);
        if (status[i] > 0)
        {
            global::tutor.bug("Attempt to start timer \"" + id
                              + "\" which was started already before.");
        }
        nStarts[i] += 1.0;
        status[i] = 1;
        tStart[i] = current;
        tStop[i] = current;
    }

    return;
}

void Timer::stop(const String& id)
{
    // Take time immediately, check if timer exists later.
    auto current = hrc::now();

    auto pos = std::find(ids.begin(), ids.end(), id);

    if (pos == ids.end())
    {
        global::tutor.bug("Attempt to stop new Timer \"" + id + "\" which was never started.");
    }
    else
    {
        const Size i = std::distance(ids.begin(), pos);
        // If timer was not started before it cannot be stopped.
        if (status[i] == 0)
        {
            global::tutor.bug("Attempt to stop timer \"" + id
                              + "\" which is currently not running.");
        }
        // Store current stop time.
        tStop[i] = current;
        status[i] = 0;
        // Accumulate total time of this timer.
        const Real dt = std::chrono::duration<Real>(tStop[i] - tStart[i]).count();
        tSum[i] += dt;
        tMin[i] = std::min(tMin[i], dt);
        tMax[i] = std::max(tMax[i], dt);
    }

    return;
}

void Timer::writeToScreen()
{
    Real total = 0.0;
    for (const Real& t : tSum) total += t;
    String sTotal = str("TIMER (TOTAL = %10.2E)", total);
    Size   l = std::max(sTotal.length(), string_tools::maxLength(ids));
    String idfmt = str("%%-%zus", l);
    String d(10, '-');
    std::cout << "\n";
    std::cout << str((idfmt + " %10s %10s %10s %10s %10s %10s\n").c_str(),
                     sTotal.c_str(),
                     "N    ",
                     "MIN   ",
                     "MAX   ",
                     "AVE   ",
                     "SUM   ",
                     "%TOTAL  ");
    std::cout << str((idfmt + "|%10s|%10s|%10s|%10s|%10s|%10s|\n").c_str(),
                     String(l, '-').c_str(),
                     d.c_str(),
                     d.c_str(),
                     d.c_str(),
                     d.c_str(),
                     d.c_str(),
                     d.c_str());
    for (Size i = 0; i < ids.size(); ++i)
    {
        std::cout << str((idfmt + " %10.0f %10.2E %10.2E %10.2E %10.2E %10.2f\n").c_str(),
                         ids[i].c_str(),
                         nStarts[i],
                         tMin[i],
                         tMax[i],
                         tSum[i] / nStarts[i],
                         tSum[i],
                         tSum[i] / total * 100.0);
    }
    std::cout << "\n";

    return;
}

namespace global
{
Timer timer;
}
