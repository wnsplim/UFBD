#ifndef Serialize_hpp
#define Serialize_hpp

#include <deque>
#include <istream>
#include <ostream>
#include <vector>

namespace Serialize {

inline void writeVec(std::ostream& os, const std::vector<double>& v) {
    os << v.size();
    for(double x : v) os << ' ' << x;
    os << '\n';
}
inline void readVec(std::istream& is, std::vector<double>& v) {
    size_t n; is >> n; v.resize(n);
    for(size_t i = 0; i < n; i++) is >> v[i];
}

inline void writeIVec(std::ostream& os, const std::vector<int>& v) {
    os << v.size();
    for(int x : v) os << ' ' << x;
    os << '\n';
}
inline void readIVec(std::istream& is, std::vector<int>& v) {
    size_t n; is >> n; v.resize(n);
    for(size_t i = 0; i < n; i++) is >> v[i];
}

inline void writeLVec(std::ostream& os, const std::vector<long>& v) {
    os << v.size();
    for(long x : v) os << ' ' << x;
    os << '\n';
}
inline void readLVec(std::istream& is, std::vector<long>& v) {
    size_t n; is >> n; v.resize(n);
    for(size_t i = 0; i < n; i++) is >> v[i];
}

inline void write2D(std::ostream& os, const std::vector<std::vector<double>>& m) {
    os << m.size() << '\n';
    for(const std::vector<double>& r : m) writeVec(os, r);
}
inline void read2D(std::istream& is, std::vector<std::vector<double>>& m) {
    size_t n; is >> n; m.resize(n);
    for(size_t i = 0; i < n; i++) readVec(is, m[i]);
}

inline void writeBoolDeque(std::ostream& os, const std::deque<bool>& d) {
    os << d.size();
    for(bool b : d) os << ' ' << (b ? 1 : 0);
    os << '\n';
}
inline void readBoolDeque(std::istream& is, std::deque<bool>& d) {
    size_t n; is >> n; d.clear();
    for(size_t i = 0; i < n; i++){ int b; is >> b; d.push_back(b != 0); }
}

}

#endif
