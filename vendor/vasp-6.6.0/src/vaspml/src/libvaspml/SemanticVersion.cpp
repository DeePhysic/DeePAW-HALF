#include "SemanticVersion.hpp"

#include "text.hpp"

#include <algorithm>  // for min
#include <stdexcept>  // for runtime_error

using namespace vaspml;

SemanticVersion::SemanticVersion(Int major, Int minor, Int patch, String suffix)
{
    update(major, minor, patch, suffix);
}

SemanticVersion::SemanticVersion(String version)
{
    update(version);
}

SemanticVersion& SemanticVersion::update(Int major, Int minor, Int patch, String suffix)
{
    major_ = major;
    minor_ = minor;
    patch_ = patch;
    suffix_ = suffix;

    processSuffix();
    checkValidity();

    return *this;
}

SemanticVersion& SemanticVersion::update(String version)
{
    parseVersionString(version, major_, minor_, patch_, suffix_);
    update(major_, minor_, patch_, suffix_);

    return *this;
}

String SemanticVersion::toString() const
{
    return std::to_string(major_) + '.' + std::to_string(minor_) + '.' + std::to_string(patch_)
         + suffix_;
}

Int SemanticVersion::major() const
{
    return major_;
}

Int SemanticVersion::minor() const
{
    return minor_;
}

Int SemanticVersion::patch() const
{
    return patch_;
}

String SemanticVersion::suffix() const
{
    return suffix_;
}

String SemanticVersion::prerelease() const
{
    return prerelease_;
}

String SemanticVersion::build() const
{
    return build_;
}

String SemanticVersion::kind() const
{
    return kind_;
}

bool SemanticVersion::operator<(const SemanticVersion& rhs) const
{
    // clang-format off
    if      (major_ < rhs.major_) return true;
    else if (major_ > rhs.major_) return false;
    if      (minor_ < rhs.minor_) return true;
    else if (minor_ > rhs.minor_) return false;
    if      (patch_ < rhs.patch_) return true;
    else if (patch_ > rhs.patch_) return false;
    // clang-format on
    // An empty pre-release string indicates a higher version.
    if (prerelease_.empty() != rhs.prerelease_.empty())
    {
        if (prerelease_.empty()) return false;
        else return true;
    }
    // Otherwise, compare lexicographically.
    if (prerelease_ < rhs.prerelease_) return true;
    else if (prerelease_ > rhs.prerelease_) return false;

    return false;
}

bool SemanticVersion::operator==(const SemanticVersion& rhs) const
{
    if (major_ != rhs.major_) return false;
    if (minor_ != rhs.minor_) return false;
    if (patch_ != rhs.patch_) return false;
    if (prerelease_ != rhs.prerelease_) return false;

    return true;
}

bool SemanticVersion::operator!=(SemanticVersion const& rhs) const
{
    return !((*this) == rhs);
}
bool SemanticVersion::operator>(SemanticVersion const& rhs) const
{
    return rhs < (*this);
}
bool SemanticVersion::operator<=(SemanticVersion const& rhs) const
{
    return !((*this) > rhs);
}
bool SemanticVersion::operator>=(SemanticVersion const& rhs) const
{
    return !((*this) < rhs);
}

void SemanticVersion::parseVersionString(String  version,
                                         Int&    major,
                                         Int&    minor,
                                         Int&    patch,
                                         String& suffix)
{
    // Initialize to defaults.
    major = 0;
    minor = 0;
    patch = 0;
    suffix = "";

    if (version.empty())
    {
        throw std::runtime_error("ERROR: Empty version string passed to SemanticVersion.");
    }

    // Major version number.
    //**********************
    Size i = version.find('.');
    // Single number version, e.g. just "2", will convert to "2.0.0".
    if (i == String::npos)
    {
        major = text::convert<Int>(version);
        return;
    }
    // Only "." without a preceding number is invalid.
    else if (i < 1)
    {
        throw std::runtime_error("ERROR: Invalid version string (early major/minor separator)");
    }
    // Read major version number.
    else major = text::convert<Int>(version.substr(0, i));
    // If string ends after ".", finish parsing.
    if (i + 1 == version.size()) return;
    // Cut away already processed version string.
    version = version.substr(i + 1);

    // Minor version number.
    //**********************
    i = version.find('.');
    // Two-number version, e.g. just "2.1", will convert to "2.1.0".
    if (i == String::npos)
    {
        minor = text::convert<Int>(version);
        return;
    }
    // Only "." without a preceding number is invalid.
    else if (i < 1)
    {
        throw std::runtime_error("ERROR: Invalid version string (early minor/patch separator)");
    }
    // Read minor version number.
    else minor = text::convert<Int>(version.substr(0, i));
    // If string ends after ".", finish parsing.
    if (i + 1 == version.size()) return;
    // Cut away already processed version string.
    version = version.substr(i + 1);

    // Patch version number.
    //**********************
    i = std::min(version.find('-'), version.find('+'));
    // No suffix, e.g. just "2.1.3".
    if (i == String::npos)
    {
        patch = text::convert<Int>(version);
        return;
    }
    // Read patch version number.
    else patch = text::convert<Int>(version.substr(0, i));
    // Cut away already processed version string and store as suffix.
    suffix = version.substr(i);

    // Split suffix and use only first part if there are spaces.
    i = suffix.find(' ');
    if (i != String::npos) { suffix = suffix.substr(0, i); }

    return;
}

void SemanticVersion::processSuffix()
{
    // Initialize substrings.
    prerelease_ = "";
    build_ = "";

    if (suffix_.empty())
    {
        kind_ = "release";
        return;
    }

    // Split suffix and use only first part if there are spaces.
    Size i = suffix_.find(' ');
    if (i != String::npos) { suffix_ = suffix_.substr(0, i); }

    if (suffix_.at(0) == '-')
    {
        i = suffix_.find('+');
        if (i != String::npos)
        {
            if (i > 1) prerelease_ = suffix_.substr(1, i - 1);
            if (i < suffix_.size() - 1) build_ = suffix_.substr(i + 1);
            kind_ = "pre-release-dev";
        }
        else
        {
            prerelease_ = suffix_.substr(1);
            kind_ = "pre-release";
        }
    }
    else if (suffix_.at(0) == '+')
    {
        build_ = suffix_.substr(1);
        kind_ = "dev";
    }
    else
    {
        throw std::runtime_error("ERROR: Invalid semantic version, "
                                 "suffix must start with \"-\" or \"+\"");
    }

    return;
}

void SemanticVersion::checkValidity() const
{
    // Check if all version numbers are greater than zero.
    String wrongVersion;
    if (major_ < 0) wrongVersion = "major";
    if (minor_ < 0) wrongVersion = "minor";
    if (patch_ < 0) wrongVersion = "patch";
    if (!wrongVersion.empty())
    {
        throw std::runtime_error("ERROR: Invalid semantic version number \"" + wrongVersion
                                 + "\", contains negative value.");
    }

    if (kind_.find("pre-release") != String::npos)
    {
        if (prerelease_.empty())
        {
            throw std::runtime_error("ERROR: Invalid semantic version suffix \"" + suffix_
                                     + "\", empty pre-release string.");
        }
    }

    if (kind_.find("dev") != String::npos)
    {
        if (build_.empty())
        {
            throw std::runtime_error("ERROR: Invalid semantic version suffix \"" + suffix_
                                     + "\", empty build metadata string.");
        }
        if (build_.find("+") != String::npos)
        {
            throw std::runtime_error("ERROR: Invalid semantic version suffix \"" + suffix_
                                     + "\", build metadata contains another \"+\" character.");
        }
    }

    return;
}
