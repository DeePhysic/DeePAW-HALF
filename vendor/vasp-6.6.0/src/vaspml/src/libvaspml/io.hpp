#ifndef IO_HPP
#define IO_HPP

#include "types.hpp"

#include <fstream>

namespace vaspml
{
struct Record;
}

namespace vaspml::io
{

/// Describes state of an ML tag (default value, user-defined,...)
enum class TagState : Int
{
    /// Tag is not set yet.
    UNSET = 1,
    /// Tag is set to default value.
    DEFAULT = 2,
    /// Tag is set to user-defined value from INCAR file.
    INCAR = 3,
    /// Tag is set to user-defined value from INCAR file (alternative tag name).
    INCAR_ALT = 4,
    /// Tag was set by user in INCAR but overridden by program, e.g., taken from ML_FF file.
    OVERRIDE = 5
};

void   open(std::fstream& strm, const String& fileName, const String& mode);
void   open(std::ifstream& strm, const String& fileName, const String& mode = "r");
void   open(std::ofstream& strm, const String& fileName, const String& mode = "w");
void   close(std::fstream& strm);
void   close(std::ifstream& strm);
void   close(std::ofstream& strm);
bool   createDirectory(String directoryPath);
Int    removeRecursively(String fileOrDirectory);
void   writeBuffer(const Buffer& buffer, std::ostream& strm);
void   readBuffer(std::istream& strm, Buffer& buffer);
String toString(std::istream& strm);
void   readMlab(Record& record, std::istream& input);
void   readMlab(Record& record, String fileName = "ML_AB");
void   writeMlab(const Record& record, std::ostream& output);
void   writeMlab(const Record& record, String fileName = "ML_AB");
void   readXyz(Record& record, std::istream& input);
void   writeXyz(const Record& record, std::ostream& output);
void   updateMlffHeaderDate(Record& rec);
void   processMlff(Record& rec, std::fstream& strm, bool doWrite);
void   readMlffAndConvertUnits(Record& record, String fileName = "ML_FF");
void   readIncar(Record& record, std::istream& strm);
void   setupFromIncar(const Record& incar, Record& setup);
void   writeMllogfile(const Record& setup, std::ostream& output);

} // namespace vaspml::io

#endif
