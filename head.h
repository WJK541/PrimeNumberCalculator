#pragma once
#include <atomic>
#include <bitset>
#include <chrono>
#include <cmath>
#include <exception>
#include <iostream>
#include <map>
#include <stdio.h>
#include <string>
#include <thread>
#include <vector>

#define DEBUG 0
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
const ull th_count = std::thread::hardware_concurrency();

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

template <typename _T>
_T gcd(_T m, _T n)
{
    while (n != 0) {
        _T t = m % n;
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
    ull length;
};
Bitvector::Bitvector(ull n, bool init)
{
#if DEBUG
    if (n == 0)
        throw runtime_error("n must be positive integer!");
#endif

    length = ((n - 1) >> 3) + 1;
    list.reserve(length);
    if (init == true) {
        for (ull i = 0; i < length; i++)
            list[i] = (uc)0xff;
    } else {
        for (ull i = 0; i < length; i++)
            list[i] = (uc)0x00;
    }
}
bool Bitvector::at(ull n)
{
    ull shang = n >> 3;
#if DEBUG
    if (shang >= length)
        throw runtime_error("Bitlist64: out of index");
#endif
    ull yushu = n & 7;
    return (list[shang] & (1 << yushu));
}
void Bitvector::setFalse(ull n)
{
    ull shang = n >> 3;
#if DEBUG
    if (shang >= length)
        throw runtime_error("Bitlist64: out of index");
#endif
    ull yushu = n & 7;
    list[shang] = list[shang] & (~(1 << yushu));
}
void Bitvector::setTrue(ull n)
{
    ull shang = n >> 3;
#if DEBUG
    if (shang >= length)
        throw runtime_error("Bitlist64: out of index");
#endif
    ull yushu = n & 7;
    list[shang] = list[shang] | (1 << yushu);
}
void Bitvector::resetTrue()
{
    for (ull i = 0; i < length; i++)
        list[i] = (uc)0xff;
}
void Bitvector::resetFalse()
{
    for (ull i = 0; i < length; i++)
        list[i] = (uc)0x00;
}
class Range
{
public:
    atomic<ull> var;
    atomic<bool> flag = true;
    ull up;
    Range(ull start, ull end)
    {
        if (start > end)
            throw runtime_error("START can not be greater than END!");
        ull r = start % 6;
        if (r == 0 or r == 1) {
            start = ((start / 6) * 6) + 1;
            flag = false;
        } else {
            flag = true;
            start = ((start / 6) * 6) + 5;
        }
        var.store(start);
        up = end;
    }
    ull get()
    {
        if (flag) {
            flag = false;
            return var.fetch_add(2);
        }
        flag = true;
        return var.fetch_add(4);
    }
};

inline ull qpow(ull base, ull exp, ull mod)
{
    ull res = 1;
    while (exp) {
        if (exp & 1)
            res = qmul(res, base, mod);
        base = qmul(base, base, mod);
        exp >>= 1;
    }
    return res;
}
inline ll qpow(ll base, ll exp, ll mod)
{
    ll res = 1;
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

bool TrialAndError(ull n, const ull *P = prime)
{
    for (ull p : prime) {
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
    for (ull p : millerrabin_prime) {
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

void PrimeListMT_Impl(ull end, vector<ull> *primes, Range *range)
{
    ull num;
    for (;;) {
        num = range->get();
        if (num > end)
            return;
        for (ull i = 2; i < 25; i++) {
            if (num % prime[i] == 0) {
                if (num == prime[i])
                    primes->push_back(num);
                goto NEXT_LOOP;
            }
        }
        if (MillerRabin(num))
            primes->push_back(num);
    NEXT_LOOP:;
    }
}
vector<ull> PrimeListMT(ull start, ull end)
{ // PrimeList_Multithreaded
    vector<ull> primes;
    if (start <= 2 and 2 <= end) {
        primes.push_back(2);
    }
    if (start <= 3 and 3 <= end) {
        primes.push_back(3);
    }
    start = max(3, start);
    Range range(start, end);
    vector<thread> threads;
    for (ull i = 0; i < 2; i++)
        threads.push_back(thread(PrimeListMT_Impl, end, &primes, &range));
    for (ull i = 0; i < 2; i++)
        threads[i].join();
    return primes;
}

void EratosthenesSieve_Impl(ull S, vector<ull> *primes, vector<vector<ull>> *result, atomic<ull> *it, ull n, ull count)
{
    ull i, j, k, l, start, p;
    Bitvector is_prime(S >> 1, true);
    for (k = ++(*it); k * S <= n; k = ++(*it)) {
        start = k * S;
        is_prime.resetTrue();
        for (l = 0; l < count; l++) {
            p = (*primes)[l];
            for (j = max(((start + p - 1) / p) | 1, p) * p - start; j < S; j += p << 1)
                is_prime.setFalse(j >> 1);
        }
        for (i = 1; i < S and start + i <= n; i += 2) {
            if (is_prime.at(i >> 1)) {
                (*result)[k].push_back(i);
            }
        }
    }
}
vector<ull> EratosthenesSieve(ull n)
{ // multithreaded blocked Eratosthenes sieve
    const ull S = 524288;
    atomic<ull> it(0);
    vector<ull> primes;
    vector<vector<ull>> result;
    for (ull i = 0; i < (n / S + 1); i++) {
        result.push_back(vector<ull>{});
    }
    if (n >= 2)
        result[0].push_back(2);
    ull nsqrt = (ull)sqrt(n), count = 0, i, j, p;
    Bitvector is_prime0((nsqrt >> 1) + 1, true);
    for (i = 3; i <= nsqrt; i += 2) {
        if (is_prime0.at(i >> 1)) {
            primes.push_back(i);
            count++;
            for (j = i * i; j <= nsqrt; j += i << 1)
                is_prime0.setFalse(j >> 1);
        }
    }
    Bitvector is_prime(S >> 1, true);
    for (i = 0; i < count; i++) {
        p = primes[i];
        for (j = p * p; j < S; j += p << 1)
            is_prime.setFalse(j >> 1);
    }
    for (i = 3; i < S and i <= n; i += 2) {
        if (is_prime.at(i >> 1)) {
            result[0].push_back(i);
        }
    }
    vector<thread> threads;
    for (i = 0; i < 2; i++)
        threads.push_back(thread(EratosthenesSieve_Impl, S, &primes, &result, &it, n, count));

    for (i = 0; i < 2; i++)
        threads[i].join();

    vector<ull> a;
    for (const auto b : result) {
        a.insert(a.end(), b.begin(), b.end());
    }
    return a;
}