#pragma once

#include <iostream>
#include <string>
#include <cassert>
#include <optional>
#include <vector>

namespace utils {

using namespace std::string_literals;

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

inline std::string subtituteStr(const std::string& text,
                                const std::string& pattern,
                                const std::string& pattern_new)
{
    typedef std::string::size_type ssize_type;
    const ssize_type n_old = pattern.length();
    const ssize_type n_new = pattern_new.length();
    ssize_type pos_begin = 0;
    ssize_type pos_found = text.find(pattern, pos_begin);
    if (pos_found == std::string::npos) {
        return text;
    }
    std::string text_new{ text };
    do {
        text_new = text_new.substr(0, pos_found) + pattern_new +
                   text_new.substr(pos_found + n_old);
        pos_begin = pos_found + n_new;
        pos_found = text_new.find(pattern, pos_begin);
    } while (pos_found != std::string::npos);
    return text_new;
}

inline void printInfo(const std::string& msg, const std::string& fname = "")
{
    bool split_newlines = (msg.find('\n') != std::string::npos);
    std::string prefix{ "[Info]: " };
    std::cout << Line_Str << '\n';
    if (not fname.empty()) {
        prefix = "[Info] (" + fname + "): ";
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
        prefix = "[Error] (" + fname + "): ";
    }
    std::cerr << prefix;
    std::cerr << (split_newlines ? subtituteStr(msg, "\n", prefix) : msg)
              << '\n';
    std::cerr << Line_Str << '\n';
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
    if (str.empty()) return false;
    char* p;
    std::strtod(trimRight(str).c_str(), &p);
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
            m_val = std::strtod(trimRight(str).c_str(), &p);
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

    const auto it_begin_last10 = vect.end() - n;
    const auto it_end = vect.end();

    for (auto it = vect.begin(); it != it_end - 2 * n; ++it) {
        if (std::equal(it_begin_last10, it_end, it)) {
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

}
