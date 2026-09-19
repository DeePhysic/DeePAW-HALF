#include "buffer.hpp"

using namespace vaspml::io::detail;

vaspml::ItemIndex vaspml::io::detail::readItemIndex(Buffer const& buffer, Size& position)
{
    ItemIndex marker;
    memcpy(&marker, &buffer[position], sizeof(ItemIndex));
    position += sizeof(ItemIndex);

    return marker;
}

void vaspml::io::detail::write(const ItemIndex input, Buffer& buffer)
{
    // Store current position in buffer.
    Size position = buffer.size();
    // Add space for single item.
    buffer.resize(buffer.size() + sizeof(ItemIndex));
    // Copy item to previous end of buffer.
    memcpy(&buffer[position], &input, sizeof(ItemIndex));

    return;
}

void vaspml::io::detail::read(Buffer const&        buffer,
                              bool&                output,
                              Size&                position,
                              [[maybe_unused]] Int unused)
{
    // Logicals in Fortran are usually 4 bytes long.
    int32_t logical;
    if (buffer.size() - position < sizeof(int32_t)) errorOutOfBuffer<int32_t>(VASPML_FLF);
    memcpy(&logical, &buffer[position], sizeof(int32_t));
    position += sizeof(int32_t);
    // Cast integer to bool.
    output = logical;

    return;
}

void vaspml::io::detail::read(Buffer const& buffer, String& output, Size& position, Int nchars)
{
    UInt64 bytesLeft = buffer.size() - position;
    if (bytesLeft <= 0) errorOutOfBuffer<String>(VASPML_FLF);

    const char* sbuf = reinterpret_cast<const char*>(&buffer[position]);
    // Either search for null terminator or use input string length.
    if (nchars < 0)
    {
        const char* spos = std::find(sbuf, sbuf + bytesLeft, '\0');
        if (spos == sbuf + bytesLeft) errorOutOfBuffer<String>(VASPML_FLF);
        position += (std::distance(sbuf, spos) + 1) * sizeof(char);
        // Construct string from char*, searches automatically for null terminator.
        output = String(sbuf);
        return;
    }
    else
    {
        if (nchars > static_cast<Int>(bytesLeft)) errorOutOfBuffer<String>(VASPML_FLF);
        position += nchars * sizeof(char);
        output = String(sbuf, nchars);
        return;
    }
}

void vaspml::io::detail::read(Buffer const& buffer, Vec1String& output, Size& position, Int nchars)
{
    Size n;
    // Either size is written in buffer or determined from output.
    if (nchars < 0) read(buffer, n, position);
    else n = output.size();

    output.resize(n);
    for (String& s : output) read(buffer, s, position, nchars);

    return;
}
