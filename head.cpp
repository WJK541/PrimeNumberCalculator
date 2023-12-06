#pragma once
#include <string>
#include <map>

typedef unsigned int uint;
typedef unsigned long long ull;
typedef unsigned char uc;
typedef long long ll;
typedef long double ld;

std::map<std::string, const std::string> test_set = {
    {"test1000", "111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111991"},
    {"test500", "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111371"},
    {"test400", "1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111373"},
    {"test300", "111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111827"},
    {"test200", "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111213"},
    {"test100", "1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111153"},
    {"test50", "11111111111111111111111111111111111111111111111119"},
    {"test40", "1111111111111111111111111111111111111117"},
    {"test30", "111111111111111111111111111191"},
    {"test20", "11111111111111111131"},
    {"test10", "1111111121"},
};
class Bitslist64
{
public:
    Bitslist64(ull n, bool init);
    ~Bitslist64();
    inline bool at(ull index);
    inline void setfalse(ull index);
    inline void settrue(ull index);
    void reset(bool a);

private:
    char *list;
    ull length;
};
Bitslist64::Bitslist64(ull n, bool init)
{
    if (n == 0)
        throw "n must be positive integer!";

    length = ((n - 1) >> 3) + 1;
    list = new char[length];
    if (init == true) {
        for (ull i = 0; i < length; i++)
            list[i] = (char)0xff;
    } else {
        for (ull i = 0; i < length; i++)
            list[i] = (char)0x00;
    }
}
Bitslist64::~Bitslist64()
{
    delete[] list;
}
inline bool Bitslist64::at(ull n)
{
    ull shang = n >> 3;
    ull yushu = n & 7;
    return (list[shang] & (1 << yushu));
}
inline void Bitslist64::setfalse(ull n)
{
    ull shang = n >> 3;
    ull yushu = n & 7;
    list[shang] = list[shang] & (~(1 << yushu));
}
inline void Bitslist64::settrue(ull n)
{
    ull shang = n >> 3;
    ull yushu = n & 7;
    list[shang] = list[shang] | (1 << yushu);
}
void Bitslist64::reset(bool a)
{
    if (a == true) {
        for (ull i = 0; i < length; i++)
            list[i] = (char)0xff;
    } else {
        for (ull i = 0; i < length; i++)
            list[i] = (char)0x00;
    }
}