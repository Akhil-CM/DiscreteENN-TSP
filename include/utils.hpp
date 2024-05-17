// -*- C++ Header -*-
#pragma once

#include <cstdio>
#include <iostream>
#include <string>
#include <cassert>
#include <cmath>
#include <optional>
#include <vector>
#include <algorithm>
#include <filesystem>

namespace utils {

namespace stdfs = std::filesystem;

using namespace std::string_literals;

const float pi_1_4{ static_cast<float>(std::atan2(1, 1)) };
const float pi_1_2{ static_cast<float>(std::atan2(1, 0)) };
const float pi_3_4{ static_cast<float>(std::atan2(1, -1)) };
const float pi{ static_cast<float>(std::atan2(0, -1)) };

// Floating point precision helpers
constexpr float epsilonf = 1e-2;
inline bool isEqual(float a, float b)
{
    return std::abs(a - b) < epsilonf;
}

inline float getRound2(float value)
{
    const float round1{ std::round(value * 1000.f) };
    const float round2{ std::round(round1/10.f) };
    return (round2/100.f);
}

inline float getRound3(float value)
{
    const float round1{ std::round(value * 10000.f) };
    const float round2{ std::round(round1/10.f) };
    return (round2/1000.f);
}

inline float getRoundN(float value, int places)
{
    const float round1 = std::round(value * std::pow(10.f, places+1));
    const float round2{ std::round(round1/10.f) };
    return (round2/std::pow(10.f, places));
}

// Console write related helpers
const std::string Line_Str = std::string{}.assign(30, '-');
const std::string Whitespace_Str = " \n\r\t\f\v";

inline void printNewLine(int count = 1, std::ostream& stream = std::cout)
{
    assert("[Error] (printNewLine): non-positive value given as count" &&
           (count > 0));
    const std::string newline_str =
        std::string{}.assign(static_cast<std::size_t>(count), '\n');
    stream << newline_str;
}

inline std::string subtituteStr(std::string text,
                                const std::string& pattern,
                                const std::string& pattern_new)
{
    typedef std::string::size_type ssize_type;
    const ssize_type n_old = pattern.length();
    const ssize_type n_new = pattern_new.length();
    ssize_type pos_found = text.find(pattern, 0);
    while (pos_found != std::string::npos) {
        text = text.substr(0, pos_found) + pattern_new +
                   text.substr(pos_found + n_old);
        pos_found = text.find(pattern, pos_found + n_new);
    } ;
    return text;
}

template<typename... Args>
inline std::string stringFmt(const std::string& fmt, Args... args)
{
    const int size = std::snprintf(nullptr, 0, fmt.c_str(), args...);
    if (size < 0) {
        throw std::runtime_error{ "[Error] (stringFmt): during formatting with:\n" + fmt };
    }
    auto buf = std::make_unique<char[]>(size + 1);
    std::snprintf(buf.get(), size + 1, fmt.c_str(), args...);
    std::string result{buf.get(), static_cast<std::string::size_type>(size)};
    return result;
}

template<typename... Args>
inline void printInfoFmt(const std::string& fmt, const std::string& fname, Args... args)
{
    std::string msg{ stringFmt(fmt, args...) };
    bool split_newlines = (msg.find('\n') != std::string::npos);
    std::string prefix{ "[Info]: " };
    std::cout << Line_Str << '\n';
    if (not fname.empty()) {
        prefix = "\n[Info] (" + fname + "): ";
    }
    std::cout << prefix;
    std::cout << (split_newlines ? subtituteStr(msg, "\n", "\n" + prefix) : msg)
              << '\n';
    std::cout << Line_Str << '\n';
}

template<typename... Args>
inline void printErrFmt(const std::string& fmt, const std::string& fname, Args... args)
{
    std::string msg{ stringFmt(fmt, args...) };
    bool split_newlines = (msg.find('\n') != std::string::npos);
    std::string prefix{ "[Error]: " };
    std::cerr << Line_Str << '\n';
    if (not fname.empty()) {
        prefix = "\n[Error] (" + fname + "): ";
    }
    std::cerr << prefix;
    std::cerr << (split_newlines ? subtituteStr(msg, "\n", "\n" + prefix) : msg)
              << '\n';
    std::cerr << Line_Str << '\n';
}

inline void printInfo(const std::string& msg, const std::string& fname = "")
{
    bool split_newlines = (msg.find('\n') != std::string::npos);
    std::string prefix{ "[Info]: " };
    std::cout << Line_Str << '\n';
    if (not fname.empty()) {
        prefix = "\n[Info] (" + fname + "): ";
    }
    std::cout << prefix;
    std::cout << (split_newlines ? subtituteStr(msg, "\n", "\n" + prefix) : msg)
              << '\n';
    std::cout << Line_Str << '\n';
}

inline void printErr(const std::string& msg, const std::string& fname = "")
{
    bool split_newlines = (msg.find('\n') != std::string::npos);
    std::string prefix{ "[Error]: " };
    std::cerr << Line_Str << '\n';
    if (not fname.empty()) {
        prefix = "\n[Error] (" + fname + "): ";
    }
    std::cerr << prefix;
    std::cerr << (split_newlines ? subtituteStr(msg, "\n", "\n" + prefix) : msg)
              << '\n';
    std::cerr << Line_Str << '\n';
}

template <typename F,
          typename std::enable_if<std::is_convertible_v<F, stdfs::path>>::type* =
              nullptr>
stdfs::path getCleanPath(const F& src)
{
    const stdfs::path tmp_src{ src };
    const stdfs::path lexical_src{ tmp_src.lexically_normal() };
    const stdfs::path abs_src{ stdfs::absolute(lexical_src) };
    return stdfs::weakly_canonical(abs_src);
}

// "Expected" type helpers
template <typename T, typename E> class Expected
{
public:
    Expected()
        : m_opt{ std::nullopt }
        , m_err{}
    {
    }
    Expected(const T& val, const E& err)
        : m_opt{ val }
        , m_err{ err }
    {
    }
    Expected(const std::optional<T>& opt, const E& err)
        : m_opt{ opt }
        , m_err{ err }
    {
    }

    bool has_value() const
    {
        return m_opt.has_value();
    }

    std::optional<T>& opt()
    {
        return m_opt;
    }
    const std::optional<T>& opt() const
    {
        return m_opt;
    }

    T& value()
    {
        return *m_opt;
    }
    const T& value() const
    {
        return *m_opt;
    }

    E& err()
    {
        return m_err;
    }
    const E& err() const
    {
        return m_err;
    }

private:
    std::optional<T> m_opt;
    E m_err;
};

template <typename T> using ErrorMsg = Expected<T, std::string>;
template <typename T> using ErrorBool = Expected<T, bool>;

inline std::string trimLeft(const std::string& text)
{
    if (text.empty()) {
        return text;
    }
    std::string::size_type pos = text.find_first_not_of(" \n\r\t\f\v");
    if (pos == std::string::npos) {
        return "";
    }
    return text.substr(pos);
}
inline std::string trimRight(const std::string& text)
{
    if (text.empty()) {
        return text;
    }
    std::string::size_type pos = text.find_last_not_of(" \n\r\t\f\v");
    if (pos == std::string::npos) {
        return "";
    }
    return text.substr(0, pos + 1);
}

inline bool isNumber(const std::string& str)
{
    if (str.empty())
        return false;
    const std::string& number_str{ trimRight(str) };
    char* p;
    std::strtod(number_str.c_str(), &p);
    return (*p == '\0');
}

class Str2Num
{
public:
    Str2Num(const std::string& str)
    {
        if (str.empty()) {
            m_val = -1.0;
            m_state = false;
        } else {
            char* p;
            const std::string& number_str{ trimRight(str) };
            m_val = std::strtod(number_str.c_str(), &p);
            m_state = (*p == '\0');
        }
    }
    bool has_value()
    {
        return m_state;
    }
    double value()
    {
        return m_val;
    }

private:
    double m_val;
    bool m_state{ false };
};

template <typename T> class MatchItem
{
public:
    MatchItem(T item)
        : m_item{ item }
    {
    }
    bool operator()(const T& item)
    {
        return (item == m_item);
    }
    T& item()
    {
        return m_item;
    }

    const T& item() const
    {
        return m_item;
    }

private:
    T m_item;
};

template <typename T>
bool hasRepeatingPattern(const std::vector<T>& vect, int n)
{
    if (vect.size() <
        static_cast<std::size_t>(
            2 * n)) { // If there aren't enough elements for a duplicate
        return false;
    }

    const auto it_begin_pattern = vect.end() - n;
    const auto it_end = vect.end();
    const auto it_end_check = it_end - (2 * n);

    for (auto it = vect.begin(); it != it_end_check; ++it) {
        if (std::equal(it_begin_pattern, it_end, it)) {
            return true;
        }
    }

    return false;
}

template <typename T>
bool vectContains(const T& item, const std::vector<T>& vect)
{
    for (const T& elmnt : vect) {
        if (elmnt == item) {
            return true;
        }
    }
    return false;
}

template <typename T> int vectFind(const T& item, const std::vector<T>& vect)
{
    const auto& it_begin{ vect.begin() };
    const auto& it_end{ vect.end() };
    const auto& it =
        std::find_if(it_begin, it_end, utils::MatchItem<T>{ item });
    return (it == it_end) ? -1 : std::distance(it_begin, it);
}

typedef std::uint32_t hash_type;

inline std::uint32_t calcCRC32(const std::string& str)
{
    const unsigned char *pData{ reinterpret_cast<const unsigned char*>(str.c_str()) };
    std::size_t ulByteCount{ str.length() };
    std::uint32_t d, ind;
    std::uint32_t acc = 0xFFFFFFFF;
    const std::uint32_t ulCrcRand32Lut[] =
    {
        0x00000000, 0x1DB71064, 0x3B6E20C8, 0x26D930AC,
        0x76DC4190, 0x6B6B51F4, 0x4DB26158, 0x5005713C,
        0xEDB88320, 0xF00F9344, 0xD6D6A3E8, 0xCB61B38C,
        0x9B64C2B0, 0x86D3D2D4, 0xA00AE278, 0xBDBDF21C
    };

    while ( ulByteCount > 0 )
    {
        ulByteCount--;
        d = *pData++;
        ind = (acc & 0x0F) ^ (d & 0x0F);
        acc = (acc >> 4) ^ ulCrcRand32Lut[ind];
        ind = (acc & 0x0F) ^ (d >> 4);
        acc = (acc >> 4) ^ ulCrcRand32Lut[ind];
    }

    return (acc ^ 0xFFFFFFFF);
}

inline hash_type getHash(const std::string& str)
{
    return calcCRC32(str);
    // return std::hash<std::string>{}(str);
}

} // namespace utils
