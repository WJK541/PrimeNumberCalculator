#pragma once
#include <atomic>
#include <bitset>
#include <chrono>
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <stdio.h>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#define DEBUG 1
#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define qmul(x, y, mod) (((__int128)x * y) % mod)
#define pi_n(n) (n / (log(n) - 1.1))

using namespace std;
using namespace chrono;

steady_clock::time_point start_time;
double elapsed_time;
#define START_TIMING start_time = steady_clock::now();
#define END_TIMING_S                                                                              \
    elapsed_time = duration_cast<miliseconds>(steady_clock::now() - start_time).count() / 1000.0; \
    printf("The calculation takes %.4f seconds\n", elapsed_time);
#define END_TIMING_MS                                                                              \
    elapsed_time = duration_cast<microseconds>(steady_clock::now() - start_time).count() / 1000.0; \
    printf("The calculation takes %.4f miliseconds\n", elapsed_time);

typedef unsigned int uint;
typedef unsigned long long ull;
typedef unsigned char uc;
typedef long long ll;
typedef long double ld;

const ull millerrabin_prime[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
const ull prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
const ull th_count = thread::hardware_concurrency();
const ull S = 2097152;

bool restart = false;
bool print_res = false;

map<string, bool (*)()> commands;

map<string, const string> test_set = {
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

template <typename T>
T gcd(T m, T n)
{
    while (n != 0) {
        T t = m % n;
        m = n;
        n = t;
    }
    return m;
}

string toBin(char c)
{
    bitset<sizeof(char) * 8> b(c);
    return b.to_string();
}

class Bitvector
{
public:
    Bitvector(ull, bool);
    bool at(ull);
    void setFalse(ull);
    void setTrue(ull);
    void resetTrue();
    void resetFalse();
    vector<uc> list;
};
Bitvector::Bitvector(ull n, bool init)
{
#if DEBUG
    if (n == 0)
        throw runtime_error("n must be positive integer!");
#endif

    ull length = ((n - 1) >> 3) + 1;
    if (init) {
        list.resize(length, (uc)0xff);
    } else {
        list.resize(length, (uc)0x00);
    }
}
bool Bitvector::at(ull n)
{
    ull shang = n >> 3;
    ull yushu = n & 7;
#if DEBUG
    int item = 0;
    try {
        item = list.at(shang) & (1 << yushu);
    } catch (exception e) {
        cout << "at: " << shang << '\n';
    }
    return item;
#else
    return (list[shang] & (1 << yushu));
#endif
}
void Bitvector::setFalse(ull n)
{
    ull shang = n >> 3;
    ull yushu = n & 7;
#if DEBUG
    int item = 0;
    try {
        item = list.at(shang) & (~(1 << yushu));
        list.at(shang) = item;
    } catch (exception e) {
        cout << "setFalse: " << shang << '\n';
    }
#else
    list[shang] = list[shang] & (~(1 << yushu));
#endif
}
void Bitvector::setTrue(ull n)
{
    ull shang = n >> 3;
    ull yushu = n & 7;
#if DEBUG
    list.at(shang) = list.at(shang) | (1 << yushu);
#else
    list[shang] = list[shang] | (1 << yushu);
#endif
}
void Bitvector::resetTrue()
{
    for (ull i = 0; i < list.size(); i++)
        list[i] = (uc)0xff;
}
void Bitvector::resetFalse()
{
    for (ull i = 0; i < list.size(); i++)
        list[i] = (uc)0x00;
}

class Range
{ // 性能极差
public:
    vector<vector<ull> *> vv;
    Range(ull start, ull end)
    {
        ull r = start % 6;
        bool flag;
        if (r == 0 or r == 1) {
            flag = false;
            start = ((start / 6) * 6) + 1;
        } else {
            flag = true;
            start = ((start / 6) * 6) + 5;
        }
        ull i = start;
        for (;;) {
            vv.push_back(new vector<ull>());
            do {
                if (flag) {
                    flag = false;
                    vv.back()->push_back(i);
                    i += 2;
                } else {
                    flag = true;
                    vv.back()->push_back(i);
                    i += 4;
                }
                if (i > end) {
                    return;
                }
            } while (i < (start + S));
            start = i;
            // delete v;
        }
    }
    ~Range()
    {
        for (auto &v1 : vv) {
            delete v1;
        }
    }
    vector<ull> *get(ull index)
    {
        if (index < vv.size())
            return vv.at(index);
        return NULL;
    }
};

template <typename _T>
inline _T qpow(_T base, _T exp, _T mod)
{
    _T res = 1;
    while (exp) {
        if (exp & 1)
            res = qmul(res, base, mod);
        base = qmul(base, base, mod);
        exp >>= 1;
    }
    return res;
}
inline ll absmod(ll x, ll mod)
{
    x %= mod;
    if (x > 0)
        return x;
    return x + mod;
}

bool TrialAndError(ull n)
{
    for (const ull &p : prime) {
        if (n % p == 0)
            return n == p;
    }
    return true;
}

bool MillerRabin(ull n)
{
    ull nm1 = n - 1; // n minus 1
    ull t = 0, s = nm1;
    while ((s & 1) == 0) {
        s >>= 1;
        t++;
    }
    // s * 2 ^ t = nm1
    ull v;
    for (const ull &p : millerrabin_prime) {
        v = qpow(p, s, n);
        if (v == 0 or v == 1 or v == nm1)
            continue;
        for (ull j = 1; j < t; j++) {
            v = qmul(v, v, n);
            if (v == nm1) {
                // v = 1;
                // break;
                goto NEXT_LOOP;
            }
            if (v == 1)
                return false;
        }
        if (v != 1)
            return false;
    NEXT_LOOP:;
    }
    return true;
}

map<ull, ull> PrimeFactorization(ull n)
{
    ull root = (ull)sqrt(n), exp = 0, i = 5;
    map<ull, ull> result = {};

#define factorize(i)     \
    while (n % i == 0) { \
        n /= i;          \
        exp++;           \
    }                    \
    if (exp) {           \
        result[i] = exp; \
        exp = 0;         \
    }                    \
    if (n == 1)          \
        return result;

    factorize(2);
    factorize(3);
    while (i <= root) {
        factorize(i);
        i += 2;
        factorize(i);
        i += 4;
    }
    if (n > 1) {
        result[n] = 1;
    }
    return result;
}

vector<ull> PrimeList64(ull start, ull end)
{
    ull i, j, root = (ull)sqrt(end) + 1;
    Bitvector is_prime((end + 1) / 2, true);
    vector<ull> primes;
    primes.reserve(pi_n(end) - pi_n(start));
    if (start <= 2 and 2 <= end)
        primes.push_back(2);
    start = max(3, start);
    for (i = 3; i < start; i += 2) {
        if (is_prime.at(i >> 1)) {
            for (j = i * i; j <= end; j += i << 1) {
                is_prime.setFalse(j >> 1);
            }
        }
    }
    for (i = start; i < root; i += 2) {
        if (is_prime.at(i >> 1)) {
            primes.push_back(i);
            for (j = i * i; j <= end; j += i << 1) {
                is_prime.setFalse(j >> 1);
            }
        }
    }
    for (i = max(root, start) | 1; i <= end; i += 2) {
        if (is_prime.at(i >> 1)) {
            primes.push_back(i);
        }
    }
    return primes;
}

void PrimeListMT_Impl(ull end, atomic<ull> *index, vector<vector<ull> *> *primes, Range *range)
{
    ull num, k;
    vector<ull> *v;
    for (;;) {
        k = index->fetch_add(1);
        v = range->get(k);
        if (v == NULL) {
            return;
        }
        for (const ull &num : *v) {
            for (ull i = 2; i < 25; i++) {
                if (num % prime[i] == 0) {
                    if (num == prime[i])
                        primes->at(k)->push_back(num);
                    goto NEXT_LOOP;
                }
            }
            if (MillerRabin(num))
                primes->at(k)->push_back(num);
        NEXT_LOOP:;
        }
    }
}
vector<vector<ull> *> PrimeListMT(ull start, ull end)
{ // PrimeList_Multithreaded
    atomic<ull> index = 0;
    vector<vector<ull> *> primes;
    for (ull i = 0; i < ((end - start) / S + 1); i++)
        primes.push_back(new vector<ull>());
    if (start <= 2 and 2 <= end)
        primes[0]->push_back(2);
    if (start <= 3 and 3 <= end)
        primes[0]->push_back(3);
    start = max(3, start);
    Range range(start, end);
    vector<thread> threads;
    for (ull i = 0; i < 32; i++)
        threads.push_back(thread(PrimeListMT_Impl, end, &index, &primes, &range));
    for (auto &th : threads) {
        th.join();
    }
    return primes;
}

void EratosthenesSieve_Impl(vector<ull> *primes, vector<vector<ull> *> *result, atomic<ull> *it, ull n)
{
    ull i, j, k, l, start, p;
    Bitvector is_prime(S >> 1, true);
    for (k = ++(*it); k * S <= n; k = ++(*it)) {
        start = k * S;
        is_prime.resetTrue();
        for (l = 0; l < primes->size(); l++) {
            p = (*primes)[l];
            for (j = max(((start + p - 1) / p) | 1, p) * p - start; j < S; j += p << 1)
                is_prime.setFalse(j >> 1);
        }
        for (i = 1; i < S and i <= n - start; i += 2) {
            if (is_prime.at(i >> 1)) {
                result->at(k)->push_back(i);
            }
        }
    }
}
vector<vector<ull> *> EratosthenesSieve(ull n)
{ // multithreaded blocked Eratosthenes sieve
    atomic<ull> it(0);
    vector<ull> primes;
    vector<vector<ull> *> result;
    for (ull i = 0; i < (n / S + 1); i++) {
        result.push_back(new vector<ull>());
        // pi_n((i + 1) * S) - pi_n(i * S)
    }
    if (n >= 2) {
        result[0]->push_back(2);
    }
    ull nsqrt = (ull)sqrt(n), i, j, p;
    Bitvector is_prime0((nsqrt >> 1) + 1, true);
    for (i = 3; i <= nsqrt; i += 2) {
        if (is_prime0.at(i >> 1)) {
            primes.push_back(i);
            for (j = i * i; j <= nsqrt; j += i << 1)
                is_prime0.setFalse(j >> 1);
        }
    }
    Bitvector is_prime(S >> 1, true);
    for (i = 0; i < primes.size(); i++) {
        p = primes[i];
        for (j = p * p; j < S; j += p << 1)
            is_prime.setFalse(j >> 1);
    }
    for (i = 3; i < S and i <= n; i += 2) {
        if (is_prime.at(i >> 1)) {
            result[0]->push_back(i);
        }
    }
    vector<thread> threads;
    const ull thread_count = min(n / S, th_count);
    // const ull thread_count = 32;
    for (i = 0; i < thread_count; i++)
        threads.push_back(thread(EratosthenesSieve_Impl, &primes, &result, &it, n));

    for (auto &th : threads)
        th.join();

    return result;
}

int _JacobiSymbolImpl(ll numerator, ll denominator)
{
    if (denominator == 1)
        return 1; //"Following the normal convention for the empty product, (Z / 1) = 1 for all `Z` "
    if (numerator == 0)
        return 0;
    numerator = numerator % denominator;
    uint countFactorsOf2 = 0;
    for (; numerator % 2 == 0 && numerator > 0; numerator /= 2)
        countFactorsOf2++;
    if (gcd(numerator, denominator) != 1)
        return 0;
    int ret = 1;
    if (countFactorsOf2 & 1) {
        uint m = denominator % 8;
        ret = ((m == 1) or (m == 7)) ? 1 : -1;
    }
    if (denominator % 4 == 3 and numerator % 4 == 3)
        ret = -ret;
    return ret * _JacobiSymbolImpl(denominator, numerator);
}
int JacobiSymbol(ll upperArgument, ll lowerArgument)
{
    if (lowerArgument % 2 == 0 or lowerArgument < 0)
        throw logic_error("lowerArgument of function `JacobiSymbol` must be a positive odd number.");
    ll denominator = lowerArgument;
    ll numerator = upperArgument;
    if (numerator < 0) {
        numerator = numerator % denominator;
        numerator += denominator;
    }
    return _JacobiSymbolImpl(numerator, denominator);
}
bool lucassequence(ll num, ll D)
{
    ll P, Q, Q_k, Q_2k, V_k, V_2k, U_k, U_2k, n;

    P = 1;
    Q = ((1 - D) / 4) % num;
    Q_k = Q;
    Q_2k = 0;
    V_k = P;
    V_2k = 2;
    U_k = 1;
    U_2k = 0;
    n = num + 1;
    vector<char> bit;
    while (n > 0) {
        if (n & 1) {
            bit.push_back('1');
        } else {
            bit.push_back('0');
        }
        n >>= 1;
    }
    /*
    for (auto j = bit.end() - 1; j >= bit.begin(); j--) {
        cout << *j;
    }
    cout << "\nQ=" << Q << ", D="
         << D << '\n'; //*/
    ll a = 1;
    for (auto i = bit.end() - 2; i >= bit.begin(); i--) {
        Q_2k = ((__int128_t)Q_k * Q_k) % num;
        U_2k = ((__int128_t)U_k * V_k) % num;
        V_2k = ((__int128_t)V_k * V_k - 2 * Q_k) % num;
        a *= 2;
        /*
        cout << "U" << a << ": " << U_2k << "\n";
        cout << "V" << a << ": " << V_2k << "\n\n"; //*/

        if (*i == '1') { // 下标+1
            a += 1;
            Q_k = ((__int128_t)Q_2k * Q) % num;

            U_k = ((__int128_t)P * U_2k + V_2k);
            if (U_k & 1)
                U_k += num;
            U_k = (U_k / 2) % num;

            V_k = ((__int128_t)D * U_2k + P * V_2k);
            if (V_k & 1)
                V_k += num;
            V_k = (V_k / 2) % num;
            /*
            cout << "U" << a << ": " << U_k << "\n";
            cout << "V" << a << ": " << V_k << "\n\n"; //*/
        } else {
            Q_k = Q_2k;
            U_k = U_2k;
            V_k = V_2k;
        }
    }
    return (U_k % num == 0) && (absmod(V_k, num) == absmod(2 * Q, num));
}
ll getD(ll n)
{
    ll D = 5;
    ll s = 1;
    ll ds, g;
    for (;;) {
        ds = D * s;
        g = gcd(ds, n);
        if (g != 1 or g != -1)
            return 0;

        if (JacobiSymbol(ds, n) == -1)
            return ds;

        D += 2;
        s = -s;
    }
}
bool Lucas(ll n)
{
    for (ll p : prime) {
        if (n % p == 0)
            return n == p;
    }
    ll D = 5;
    ll s = 1;
    ll ds, g;
    for (;;) {
        ds = D * s;
        g = gcd(ds, n);
        if (abs(g) != 1)
            return 0;

        if (JacobiSymbol(ds, n) == -1)
            break;

        D += 2;
        s = -s;
    }
    return lucassequence(n, ds);
}

template <typename T>
void output(T msg)
{
    cout << msg;
}
template <>
void output<bool>(bool msg)
{
    if (msg) {
        printf("YES\n");
    } else {
        printf("NO\n");
    }
}
template <>
void output<vector<ull>>(vector<ull> msg)
{
    if (print_res) {
        for (const auto &i : msg) {
            printf("%llu ", i);
        }
    }
    cout << '\n'
         << msg.size() << '\n';
}
template <>
void output<vector<vector<ull> *>>(vector<vector<ull> *> msg)
{
    ull count = 0;
    for (auto &i : msg) {
        count += i->size();
        if (not print_res) {
            delete i;
            continue;
        }
        for (const auto &j : *i) {
            printf("%llu ", j);
        }
        delete i;
    }
    cout << '\n'
         << count << '\n';
}
