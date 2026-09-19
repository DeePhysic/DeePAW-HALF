#include "io_detail.hpp"

#include "StringAlternative.hpp"
#include "rec.hpp"

#include <algorithm> // for find, max, transform
#include <cctype>    // for toupper
#include <iterator>  // for distance
#include <regex>     // for regex_replace, regex
#include <sstream>   // for basic_ostream, operator<<, basic_is...
#include <stdexcept> // for runtime_error

using namespace vaspml;

Size io::detail::bytesLeftInStream(std::istream& strm)
{
    std::streampos pos = strm.tellg();
    strm.seekg(0, std::ios::end);
    std::streamsize len = strm.tellg() - pos;
    strm.seekg(pos);

    return len;
}

void io::detail::checkLine(const String& line, const String& search, const Int n)
{
    if (!string_tools::checkForSubstring(line, search))
    {
        throw std::runtime_error("ERROR: Invalid file, line " + std::to_string(n)
                                 + ": expected line containing string\n   \"" + search
                                 + "\"\nbut got instead:\n   \"" + line + "\"\n");
    }

    return;
}

Int io::detail::readLine(std::istream& input, String& buffer, Int skip, bool expectEOF)
{
    for (Int i = 0; i < skip + 1; ++i)
    {
        if (!getline(input, buffer))
        {
            if (expectEOF) return i;
            else throw std::runtime_error("ERROR: Unexpected end of file.");
        }
    }

    return skip + 1;
}

void io::detail::processMlffHeader(Buffer& buffer, Size& position, bool doWrite, Record& rec)
{
    bool headerPresent = false;

    if (doWrite) headerPresent = rec.contains("header");
    else
    {
        String firstWord;
        io::detail::read(buffer, firstWord, position, 5);
        if (firstWord == "ML_FF") headerPresent = true;
        position = 0;
    }

    if (!headerPresent) return;

    String           header;
    const Vec1String headerSort = {"date",
                                   "ML_LFAST",
                                   "ML_DESC_TYPE",
                                   "types",
                                   "training_structures",
                                   "local_reference_cfgs",
                                   "descriptors",
                                   "ML_IALGO_LINREG",
                                   "ML_RCUT1",
                                   "ML_RCUT2",
                                   "ML_W1",
                                   "ML_SION1",
                                   "ML_SION2",
                                   "ML_LMAX2",
                                   "ML_MRB1",
                                   "ML_MRB2",
                                   "ML_IWEIGHT",
                                   "ML_WTOTEN",
                                   "ML_WTIFOR",
                                   "ML_WTSIF"};
    if (doWrite)
    {
        Record&            hRec = rec.dget<ShRec>("header");
        String             headerInit = "ML_FF " + hRec.cget<String>("_version") + " binary ";
        std::ostringstream headerStream;
        headerStream << headerInit;
        if (!hRec.contains("_sort")) hRec.put<Vec1String>("_sort", headerSort);
        rec::toJson(hRec, headerStream, -1);
        header = headerStream.str();
        if (header.length() >= constants::mlffHeaderSize - 1)
        {
            global::tutor.warning(
                "JSON part of ML_FF header would be too large, its output will be omitted.");
            header = headerInit + "{ }";
        }
        header.resize(constants::mlffHeaderSize, ' ');
        header.back() = '\n';
    }
    io::detail::process(buffer, header, position, doWrite, constants::mlffHeaderSize);
    if (!doWrite)
    {
        rec.add("header", ItemIndex::SHREC);
        Record& hRec = rec.dget<ShRec>("header");

        Size recStart = std::distance(header.begin(), std::find(header.begin(), header.end(), '{'));

        String     headerIntroFull(header.begin(), header.begin() + recStart);
        Vec1String headerIntro = string_tools::split(headerIntroFull, " ");
        if (headerIntro.size() < 3 || headerIntro[0] != "ML_FF" || headerIntro[2] != "binary")
        {
            throw std::runtime_error("Malformatted ML_FF header intro detected: \""
                                     + headerIntroFull + "\".");
        }
        hRec.put("_version", headerIntro[1]);

        rec::fromJson(header.substr(recStart), hRec);
        hRec.put<Vec1String>("_sort", headerSort);
    }

    return;
}

void io::detail::setupKnownTags(Record& record)
{
    record.put<String>("_type", "incar");
    record.put<Vec1String>("_sort",
                           Vec1String{"tags", "duplicates", "deprecated", "tag_state_alt"});
    // String message summarizing how many tags were found/processed.
    record.put<String>("message", "");
    // Contains all INCAR tags found in file with values as strings.
    record.add("strings", ItemIndex::SHREC);
    // Subset of all INCAR tags for which types are known, values as actual types.
    record.add("tags", ItemIndex::SHREC);
    // Map of all known INCAR tags to their types (string representation of type, e.g. "Vec1Int").
    record.add("types", ItemIndex::SHREC);
    // Map of all known INCAR tags to their one-line descriptions.
    record.add("description", ItemIndex::SHREC);
    // List of all duplicated INCAR tags found (main and alternative tags are counted as the same
    // tag).
    record.add("duplicates", ItemIndex::VEC1STRING);
    // List of all alternative (deprecated) INCAR tags spellings found (even if tag with main
    // spelling) was already processed before.
    record.add("deprecated", ItemIndex::VEC1STRING);
    // List of all processed INCAR tags with alternative (deprecated) spellings found (later used
    // for "I" vs "i" tag state specifier in ML_LOGFILE).
    record.add("tag_state_alt", ItemIndex::VEC1STRING);

    // Store type assignments and tag descriptions.
    Vec2String known = constants::tagInfo;
    Record&    typesRecord = record.dget<ShRec>("types");
    Record&    descriptionRecord = record.dget<ShRec>("description");
    Size       maxDescriptionLength = 0;
    for (Vec1String& k : known)
    {
        maxDescriptionLength = std::max(maxDescriptionLength, k.back().size());
        descriptionRecord.put<String>(k.front(), k.back());
        // Remove description from the end.
        k.pop_back();
        typesRecord.put<String>(k.front(), k.back());
        // Remove type information from the end.
        k.pop_back();
    }
    record.put<Int>("max_description_length", maxDescriptionLength);

    // Store alternative tag data.
    record.add("tag_alternatives", ItemIndex::SHREC);
    record.get<ShRec>("tag_alternatives") = assignOrMakeRecord(nullptr, dataMapStringAlternative());
    StringAlternative sa(record.get<ShRec>("tag_alternatives"));
    sa.add(known);

    return;
}

void io::detail::readIncarToStrings(Record& record, const String& input)
{
    bool                   closeGroup = false;
    String                 group = "";
    String::const_iterator pos = input.begin();

    rec::checkType(record, "incar", flf(VASPML_FLF));
    Record&     stringsRecord = record.dget<ShRec>("strings");
    Vec1String& order = record.add("tag_order", ItemIndex::VEC1STRING).get<Vec1String>("tag_order");
    Vec1String& dup = record.get<Vec1String>("duplicates");
    Vec1String& dep = record.get<Vec1String>("deprecated");
    Vec1String& alt = record.get<Vec1String>("tag_state_alt");
    StringAlternative sa(record.get<ShRec>("tag_alternatives"));

    while (true)
    {
        // Store beginning of line.
        String::const_iterator old = pos;
        // Search for end of current line.
        pos = find(input, pos, "\n\r", flf(VASPML_FLF), false);
        // If at end of string, exit loop.
        if (pos == input.end()) break;
        // Store current line.
        String line(old, pos);
        // Advance main iterator to next line.
        pos++;
        // Search line for comment characters, remove comment and trim remaining string.
        String::const_iterator posComment = find(line, line.begin(), "!#", flf(VASPML_FLF), false);
        if (posComment != line.end()) line.erase(posComment, line.end());
        line = string_tools::trim(line);
        // If original line contains only comments and/or spaces the line variable is now empty.
        if (line.empty()) continue;
        // Search for group opening.
        String::const_iterator posBrace = find(line, line.begin(), "{", flf(VASPML_FLF), false);
        if (posBrace != line.end())
        {
            if (!group.empty())
            {
                global::tutor.error(
                    "Invalid INCAR: Found tag group symbol \"{\" but another group (\"" + group
                    + "\") is already open.");
            }
            group = String(line.cbegin(), posBrace - 1);
            group = string_tools::trim(group) + "/";
            std::transform(group.begin(),
                           group.end(),
                           group.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            if (group.empty())
            {
                global::tutor.error(
                    "Invalid INCAR: Found tag group symbol \"{\" but group name is empty.");
            }
            // After consuming group info, set line to trailing string.
            line.erase(line.cbegin(), posBrace + 1);
            line = string_tools::ltrim(line);
        }
        // If group is closed in this line
        if (line.back() == '}')
        {
            if (group.empty())
            {
                global::tutor.error(
                    "Invalid INCAR: Found tag group close symbol \"}\" but no group was opened.");
            }
            closeGroup = true;
            line.pop_back();
            line = string_tools::rtrim(line);
        }
        // Check if there is at least one tag-value pair in this line.
        String::const_iterator posEqual = find(line, line.begin(), "=", flf(VASPML_FLF), false);
        if (posEqual != line.end())
        {
            // Split string into multiple tag-value substrings at semicolon.
            Vec1String sub = string_tools::split(string_tools::trim(line), ";");
            Vec1String tv;
            for (Size i = 0; i < sub.size(); ++i)
            {
                String& s = sub[i];
                // Split substring into tag-value pair.
                tv = string_tools::split(s, "=");
                if (tv.size() != 2)
                {
                    global::tutor.error("Invalid INCAR: expected tag-value pair but got \"" + s
                                        + "\".");
                }
                // Check for empty tag or value string, remove whitespaces.
                for (auto& itv : tv)
                {
                    itv = string_tools::trim(itv);
                    if (itv.empty())
                    {
                        global::tutor.error("Invalid INCAR: empty tag or value string in \"" + s
                                            + "\".");
                    }
                }
                if (tv[1].front() == '"')
                {
                    // Remove opening double quotes from value.
                    tv[1].erase(tv[1].begin());
                    // Check if multi-line string actually ends in same line.
                    if (!tv[1].empty() && tv[1].back() == '"') tv[1].pop_back();
                    else
                    {
                        // Search for next occurence of double quotes in rest of file.
                        String::const_iterator endOfMulti = find(input, pos, "\"", flf(VASPML_FLF));
                        // Add the additional lines (add back previously deleted newline).
                        tv[1].append('\n' + String{pos, endOfMulti});
                        // Set current position to character after closing double quotes.
                        pos = endOfMulti + 1;
                    }
                }
                if (tv[1].back() == '\\')
                {
                    // Remove backslash.
                    tv[1].pop_back();
                    while (true)
                    {
                        // Search for end of next line.
                        String::const_iterator endOfLine =
                            find(input, pos, "\n\r", flf(VASPML_FLF));
                        // On Windows newline could be "\r\n". If so, advance one character.
                        if (*endOfLine == '\r' && endOfLine++ != input.end())
                        {
                            if (*(endOfLine + 1) == '\n') endOfLine++;
                        }
                        String                 contLine{pos, endOfLine};
                        String::const_iterator posBackSlash =
                            find(contLine, contLine.begin(), "\\", flf(VASPML_FLF), false);
                        if (posBackSlash != contLine.end())
                        {
                            tv[1].append(contLine.cbegin(), posBackSlash);
                        }
                        else
                        {
                            posComment =
                                find(contLine, contLine.begin(), "!#", flf(VASPML_FLF), false);
                            tv[1].append(contLine.cbegin(), posComment);
                            tv[1] = std::regex_replace(tv[1], std::regex("\\n"), "");
                            tv[1] = std::regex_replace(tv[1], std::regex("\\r"), "");
                            break;
                        }
                        pos = endOfLine + 1;
                    }
                }

                String tag = tv[0];
                String main = sa(tag);
                if (tag != main && std::find(dep.begin(), dep.end(), tag) == dep.end())
                {
                    dep.push_back(tag);
                }
                if (!stringsRecord.contains(group + main))
                {
                    order.push_back(group + main);
                    stringsRecord.put<String>(group + main, tv[1]);
                    if (tag != main) alt.push_back(group + main);
                }
                else
                {
                    if (std::find(dup.begin(), dup.end(), group + main) == dup.end())
                    {
                        dup.push_back(group + main);
                    }
                }
            }
        }
        if (closeGroup)
        {
            group = "";
            closeGroup = false;
        }
    }

    return;
}

void io::detail::processIncarTags(Record& record)
{
    String message;

    rec::checkType(record, "incar", flf(VASPML_FLF));
    const Record& stringsRecord = record.dget<ShRec>("strings");
    const Record& typesRecord = record.dget<ShRec>("types");
    Record&       tags = record.dget<ShRec>("tags");

    // Fill record with processed tags (key-value pairs with correct types).
    // Note: Tag could contain group prefix, e.g. "my_group/MYTAG".
    for (String gtag : stringsRecord.keys())
    {
        // Extract actual tag if inside group.
        String tag = string_tools::split(gtag, "/").back();
        // Check if type is known for this tag.
        if (typesRecord.contains(tag))
        {
            tags.add(gtag, typesRecord.cget<String>(tag));
            ItemIndex     i = tags.itemIndexOf(gtag);
            const String& s = stringsRecord.cget<String>(gtag);
            if (i == ItemIndex::REAL) tags.get<Real>(gtag) = convertIncar<Real>(s, gtag);
            else if (i == ItemIndex::INT) tags.get<Int>(gtag) = convertIncar<Int>(s, gtag);
            else if (i == ItemIndex::STRING) tags.get<String>(gtag) = convertIncar<String>(s, gtag);
            else if (i == ItemIndex::BOOL) tags.get<bool>(gtag) = convertIncar<bool>(s, gtag);
            else if (i == ItemIndex::VEC1REAL)
            {
                tags.get<Vec1Real>(gtag) = convertIncar<Vec1Real>(s, gtag);
            }
            else if (i == ItemIndex::VEC1INT)
            {
                tags.get<Vec1Int>(gtag) = convertIncar<Vec1Int>(s, gtag);
            }
            else
            {
                global::tutor.error("Reading INCAR tags of type \"" + tags.typeOf(gtag)
                                    + "\" is not implemented.");
            }
        }
    }

    message += "Identified " + std::to_string(stringsRecord.size()) + " INCAR tags.\n";
    message +=
        "Processed " + std::to_string(tags.size()) + " ML-related INCAR tags with known type.\n";
    Vec1String& dep = record.get<Vec1String>("deprecated");
    if (!dep.empty())
    {
        StringAlternative sa(record.get<ShRec>("tag_alternatives"));
        message += "WARNING: Deprecated input tags were found, please consider the replacements:\n";
        for (String d : dep) message += "* " + d + " => " + sa(d) + "\n";
    }
    Vec1String& dup = record.get<Vec1String>("duplicates");
    if (!dup.empty())
    {
        message +=
            "WARNING: Duplicated input tags were found, only the first occurence was processed:\n";
        for (String d : dup) message += "* " + d + "\n";
    }
    record.get<String>("message").append(message);

    return;
}

void io::detail::addVaspInterfaceDataToSetup(Record& record,
                                             Int     ntyp,
                                             Int     nsw,
                                             Int     ibrion,
                                             bool    mlabExists,
                                             bool    mlffExists)
{
    record.put<Int>("NTYP", ntyp);
    record.put<Int>("NSW", nsw);
    record.put<Int>("IBRION", ibrion);
    record.put<bool>("LMLABEXIST", mlabExists);
    record.put<bool>("LMLFFEXIST", mlffExists);

    return;
}

void io::detail::writeMllogfileKernel(const Record& setup, std::ostream& output)
{
    const Record& incar = setup.dcget<ShRec>("incar");

    output << mllogfileTag("ML_LFAST", setup);
    output << mllogfileTag("ML_LIB", setup);
    output << mllogfileTag("ML_PALGO", setup);
    output << mllogfileTag("ML_NCSHMEM", setup);
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Descriptor settings");
    output << "Radial descriptors:\n";
    output << "-------------------\n";
    output << mllogfileTag("ML_RCUT1", setup);
    output << mllogfileTag("ML_SION1", setup);
    output << mllogfileTag("ML_MRB1", setup);
    output << "\n";
    output << "Angular descriptors:\n";
    output << "--------------------\n";
    output << mllogfileTag("ML_DESC_TYPE", setup);
    output << mllogfileTag("ML_RCUT2", setup);
    output << mllogfileTag("ML_SION2", setup);
    output << mllogfileTag("ML_MRB2", setup);
    output << mllogfileTag("ML_LMAX2", setup);
    output << mllogfileTag("ML_LAFILT2", setup);
    output << mllogfileTag("ML_AFILT2", setup);
    output << mllogfileTag("ML_IAFILT2", setup);
    output << mllogfileTag("ML_LSPARSDES", setup);
    output << mllogfileTag("ML_NRANK_SPARSDES", setup);
    output << mllogfileTag("ML_RDES_SPARSDES", setup);
    output << mllogfileTag("ML_LSIC", setup);
    output << mllogfileTag("ML_BASIS_TYPE2", setup);
    const Int& basis_type2 = setup.cget<Int>("ML_BASIS_TYPE2");
    if (basis_type2 == 4 || basis_type2 == 5) output << mllogfileTag("ML_MRB2_MAX", setup);
    const Int& desc_type = setup.cget<Int>("ML_DESC_TYPE");
    if (desc_type == 2 || desc_type == 4 || desc_type == 7 || desc_type == 9)
    {
        output << mllogfileTag("ML_LMAX2_RED", setup);
        output << mllogfileTag("ML_MRB2_RED", setup);
        output << mllogfileTag("ML_DESC_RATIO_DUAL", setup);
        output << mllogfileTag("ML_DESC_FACTORE_TESTE", setup);
        output << mllogfileTag("ML_BASIS_TYPE2_RED", setup);
        const Int& basis_type2_red = setup.cget<Int>("ML_BASIS_TYPE2_RED");
        if (basis_type2_red == 4 || basis_type2_red == 5)
        {
            output << mllogfileTag("ML_MRB2_MAX_RED", setup);
        }
    }
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Kernel settings");
    output << mllogfileTag("ML_W1", setup);
    output << mllogfileTag("ML_NHYP", setup);
    output << mllogfileTag("ML_LSUPERVEC", setup);
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Error estimation");
    output << mllogfileTag("ML_ICRITERIA", setup);
    output << mllogfileTag("ML_IUPDATE_CRITERIA", setup);
    output << mllogfileTag("ML_CTIFOR", setup);
    output << mllogfileTag("ML_SCLC_CTIFOR", setup);
    output << mllogfileTag("ML_CSIG", setup);
    output << mllogfileTag("ML_CSLOPE", setup);
    output << mllogfileTag("ML_CX", setup);
    output << mllogfileTag("ML_CDOUB", setup);
    output << mllogfileTag("ML_CALGO", setup);
    output << mllogfileTag("ML_ESTBLOCK", setup);
    output << mllogfileTag("ML_SF_LIMIT", setup);
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Sparsification and regression");
    output << mllogfileTag("ML_EPS_LOW", setup);
    output << mllogfileTag("ML_EPS_REG", setup);
    output << mllogfileTag("ML_LBASIS_DISCARD", setup);
    output << mllogfileTag("ML_IALGO_LINREG", setup);
    output << mllogfileTag("ML_ISVD", setup);
    output << mllogfileTag("ML_IREG", setup);
    output << mllogfileTag("ML_SIGV0", setup);
    output << mllogfileTag("ML_SIGW0", setup);
    output << mllogfileTag("ML_SAVECMAT", setup);
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Weighting and sampling");
    output << mllogfileTag("ML_WTOTEN", setup);
    output << mllogfileTag("ML_WTIFOR", setup);
    output << mllogfileTag("ML_WTSIF", setup);
    output << mllogfileTag("ML_MHIS", setup);
    output << mllogfileTag("ML_NMDINT", setup);
    output << mllogfileTag("ML_IWEIGHT", setup);
    output << mllogfileTag("ML_LUSE_NAMES", setup);
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Static array sizes");
    output << mllogfileTag("ML_MB", setup);
    String tmp("Maximum number of local reference configurations in memory (max. buffer size "
               "before sparsification)");
    tmp.resize(incar.cget<Int>("max_description_length"), ' ');
    output << mllogfileTagManual(tmp, setup.cget<Int>("ML_MB"));
    output << mllogfileTag("ML_MB_MIN", setup);
    output << mllogfileTag("ML_MCONF", setup);
    output << mllogfileTag("ML_MCONF_NEW", setup);
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Special features");
    output << mllogfileTag("ML_LHEAT", setup);
    output << mllogfileTag("ML_LCOUPLE", setup);
    output << mllogfileTag("ML_NATOM_COUPLED", setup);
    output << mllogfileTag("ML_ICOUPLE", setup, "Ion index");
    output << mllogfileTag("ML_RCOUPLE", setup);
    output << mllogfileTag("ML_LEMPPOT", setup);
    if (setup.cget<bool>("ML_LEMPPOT"))
    {
        output << mllogfileTag("ML_EMPPOT_RCUT", setup);
        output << mllogfileTag("ML_SRPOT_B0", setup);
        output << mllogfileTag("ML_SRPOT_S0", setup);
        output << mllogfileTag("ML_SRPOT_N0", setup);
        output << mllogfileTag("ML_MOPOT_NM", setup);
        output << mllogfileTag("ML_MOPOT_DM", setup, "ML_MOPOT_DM");
        output << mllogfileTag("ML_MOPOT_RM", setup, "ML_MOPOT_RM");
        output << mllogfileTag("ML_MOPOT_RKM", setup, "ML_MOPOT_RKM");
        output << mllogfileTag("ML_MOPOT_IJM", setup, "ML_MOPOT_IJM");
    }
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Reference energies");
    output << mllogfileTag("ML_ISCALE_TOTEN", setup);
    Vec1String types = setup.cget<Vec1String>("TYPE");
    Size       maxChars = string_tools::maxLength(types);
    for (auto& t : types)
    {
        if (t.size() < maxChars) t.resize(maxChars, ' ');
    }
    output << mllogfileTag("ML_EATOM_REF", setup, "Element ", types);
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Input/output options");
    output << mllogfileTag("ML_LEATOM", setup);

    return;
}

void io::detail::writeMllogfileGrace(const Record& setup, std::ostream& output)
{
    const Record& incar = setup.dcget<ShRec>("incar");

    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Force field settings");
    String tmp("Cutoff radius");
    tmp.resize(incar.cget<Int>("max_description_length"), ' '),
        output << mllogfileTagManual(tmp, setup.cget<Real>("ML_RCUT1"));
    output << "\n";
    output << "\n";
    output << mllogfileSection('-', "Input/output options");
    output << mllogfileTag("ML_GRACE_MODEL", setup);

    return;
}

String io::detail::mllogfileSection(char fill, String title)
{
    String result(140, fill);
    result += '\n';

    if (fill == '-' && !title.empty()) result = title + "\n" + result;
    else if (fill == '*' && !title.empty())
    {
        title = ' ' + title + ' ';
        result = result.replace(1, title.size(), title) + '\n';
    }

    return result;
}

String io::detail::mllogfileTag(String        tag,
                                const Record& setup,
                                String        vectorCommon,
                                Vec1String    vectorIndex)
{
    const Record& incar = setup.dcget<ShRec>("incar");
    const Record& desc = incar.dcget<ShRec>("description");
    const Record& states = setup.dcget<ShRec>("states");
    const Size    maxValueLength = 13;
    const String  l = std::to_string(maxValueLength);

    String     result = desc.cget<String>(tag);
    const Size maxDescriptionLength = incar.cget<Int>("max_description_length");
    result.resize(maxDescriptionLength, ' ');
    result += " : ";
    ItemIndex i = setup.itemIndexOf(tag);
    if (i == ItemIndex::REAL) result += str(("%" + l + ".5E").c_str(), setup.cget<Real>(tag));
    else if (i == ItemIndex::INT) result += str(("%" + l + "d").c_str(), setup.cget<Int>(tag));
    else if (i == ItemIndex::STRING)
    {
        String tmp = setup.cget<String>(tag);
        if (tmp.size() < maxValueLength) tmp.insert(0, maxValueLength - tmp.size(), ' ');
        else if (tmp.size() > maxValueLength)
        {
            result += "see below...\n^^^ : ";
            if (tmp.size() < maxDescriptionLength + maxValueLength - 3)
            {
                tmp.insert(0, maxDescriptionLength + maxValueLength - 3 - tmp.size(), ' ');
            }
        }
        result += tmp;
    }
    else if (i == ItemIndex::BOOL)
    {
        String tmp;
        if (setup.cget<bool>(tag)) tmp = "T";
        else tmp = "F";
        if (tmp.size() < maxValueLength) tmp.insert(0, maxValueLength - tmp.size(), ' ');
        result += tmp;
    }
    else if (i == ItemIndex::VEC1REAL)
    {
        const Vec1Real& vec = setup.cget<Vec1Real>(tag);
        result += str(("%" + l + "zu").c_str(), vec.size());
    }
    else if (i == ItemIndex::VEC1INT)
    {
        const Vec1Int& vec = setup.cget<Vec1Int>(tag);
        result += str(("%" + l + "zu").c_str(), vec.size());
    }
    else
    {
        global::tutor.error("Writing INCAR tags of type \"" + setup.typeOf(tag)
                            + "\" is not implemented.");
    }
    result += " " + vaspml::toString(static_cast<TagState>(states.cget<Int>(tag))) + " ";

    if (i == ItemIndex::VEC1REAL || i == ItemIndex::VEC1INT) result += "SIZE( " + tag + " )";
    else
    {
        String tmp = tag;
        tmp.resize(24, ' ');
        result += tmp;
    }
    result += "\n";

    if (i == ItemIndex::VEC1REAL)
    {
        const Vec1Real& vec = setup.cget<Vec1Real>(tag);
        for (Size i = 0; i < vec.size(); ++i)
        {
            String tmp = vectorCommon;
            if (i < vectorIndex.size()) tmp += vectorIndex[i];
            if (tmp.size() < maxDescriptionLength)
            {
                tmp.insert(0, maxDescriptionLength - tmp.size(), ' ');
            }
            tmp += " : ";
            tmp += str(("%" + l + ".5E").c_str(), vec[i]);
            tmp += " " + vaspml::toString(static_cast<TagState>(states.cget<Int>(tag))) + " ";
            tmp += tag + "(" + std::to_string(i + 1) + ")\n";
            result += tmp;
        }
    }
    else if (i == ItemIndex::VEC1INT)
    {
        const Vec1Int& vec = setup.cget<Vec1Int>(tag);
        for (Size i = 0; i < vec.size(); ++i)
        {
            String tmp = vectorCommon;
            if (i < vectorIndex.size()) tmp += vectorIndex[i];
            if (tmp.size() < maxDescriptionLength)
            {
                tmp.insert(0, maxDescriptionLength - tmp.size(), ' ');
            }
            tmp += " : ";
            tmp += str("%13d", vec[i]);
            tmp += str(("%" + l + "d").c_str(), vec[i]);
            tmp += " " + vaspml::toString(static_cast<TagState>(states.cget<Int>(tag))) + " ";
            tmp += tag + "(" + std::to_string(i + 1) + ")\n";
            result += tmp;
        }
    }

    return result;
}
