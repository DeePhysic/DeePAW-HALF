#ifndef SEMANTICVERSION_HPP
#define SEMANTICVERSION_HPP

#include "types.hpp"

namespace vaspml
{

/***************************************************************************************************
 * Implementation of semantic versioning rules (version 2.0.0).
 *
 * For detailed rules description see https://semver.org/.
 * @note This implementation does not recognize dot-separated identifiers of pre-release and build
 * metadata strings as separate tokens. Hence, for comparison only the full strings are used.
 **************************************************************************************************/
class SemanticVersion
{
  public:
    SemanticVersion() = default;
    explicit SemanticVersion(Int major, Int minor = 0, Int patch = 0, String suffix = "");
    explicit SemanticVersion(String version);
    SemanticVersion& update(Int major, Int minor = 0, Int patch = 0, String suffix = "");
    SemanticVersion& update(String version);
    Int              major() const;
    Int              minor() const;
    Int              patch() const;
    String           suffix() const;
    String           prerelease() const;
    String           build() const;
    String           kind() const;
    String           toString() const;
    bool             operator<(const SemanticVersion& rhs) const;
    bool             operator==(const SemanticVersion& rhs) const;
    bool             operator!=(const SemanticVersion& rhs) const;
    bool             operator>(const SemanticVersion& rhs) const;
    bool             operator<=(const SemanticVersion& rhs) const;
    bool             operator>=(const SemanticVersion& rhs) const;
    static void      parseVersionString(String  version,
                                        Int&    major,
                                        Int&    minor,
                                        Int&    patch,
                                        String& suffix);

  private:
    Int    major_ = 0;
    Int    minor_ = 0;
    Int    patch_ = 0;
    String suffix_ = "";
    String prerelease_ = "";
    String build_ = "";
    String kind_ = "release";

    void processSuffix();
    void checkValidity() const;
};

} //namespace vaspml

#endif
