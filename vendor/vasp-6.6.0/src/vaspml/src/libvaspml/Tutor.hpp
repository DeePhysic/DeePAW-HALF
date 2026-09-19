#ifndef TUTOR_HPP
#define TUTOR_HPP

#include "types.hpp"

namespace vaspml
{

class Tutor
{
  public:
    // constructor for exception handler
    Tutor() {};

    // next are the different kind of program exceptions
    void warning(const String message) const;

    void error(const String message) const;

    void bug(const String message) const;

  private:
    // define headers and footers for handmane exception messages here
    const String headerBug = "\n"
                             "                _     ____    _    _    _____     _ \n"
                             "               | |   |  _ \\  | |  | |  / ____|   | |\n"
                             "               | |   | |_) | | |  | | | |  __    | |\n"
                             "               |_|   |  _ <  | |  | | | | |_ |   |_|\n"
                             "                _    | |_) | | |__| | | |__| |    _ \n"
                             "               (_)   |____/   \\____/   \\_____|   (_)\n"
                             "\n";

    const String footerBug = "\n"
                             "\n"
                             "If you are not a developer, you should not encounter this problem.\n"
                             "Please submit a bug report.\n"
                             "\n";

    const String headerError =
        "\n\n"
        "EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR      ###     ###     ###\n"
        "E        R     R  R     R  O     O  R     R     ###     ###     ###\n"
        "E        R     R  R     R  O     O  R     R     ###     ###     ###\n"
        "EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR       #       #       # \n"
        "E        R   R    R   R    O     O  R   R                          \n"
        "E        R    R   R    R   O     O  R    R      ###     ###     ###\n"
        "EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###\n"
        "\n";

    const String footerError = "\n"
                               "\n"
                               "  ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ... BYE!!! <----\n"
                               "\n";

    const String headerWarning = "\n\n"
                                 "      W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!\n"
                                 "      W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!\n"
                                 "      W    W  A    A  R    R  N N  N  II  N N  N  G       !!!\n"
                                 "      W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   ! \n"
                                 "      WW  WW  A    A  R   R   N   NN  II  N   NN  G    G     \n"
                                 "      W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!\n"
                                 "\n";
    const String footerWarning = "\n \n";
};

namespace global
{
inline Tutor tutor;
}

} // namespace vaspml

#endif
