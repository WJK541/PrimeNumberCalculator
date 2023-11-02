#include "BigInteger.cpp"
#include <atomic>
#include <chrono>
#include <cmath>
#include <gmpxx.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <thread>
#include <vector>

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace chrono;

constexpr auto invaild_mode = "Invalid mode!!!\n";
constexpr auto help_info = "\n[Commands]:\n/h  Get commands\n/m  Re-select mode\n/r  Turn on/off output result\n/q  Quit the program\n";
constexpr auto mode_info = "1)Generate prime number list-1    2)Primality test    3)Factorization    4)Generate prime number list-2\n5)Eratosthenes sieve    6)Miller-Rabin    7)Debug-1    8)Debug-2\n";
constexpr auto choose_mode_info = "\n[Select mode]:";
system_clock::time_point start_time;
double used_time;
bool restart = false;
bool print_res = false;

typedef unsigned int uint;
typedef unsigned long long ull;
typedef unsigned char uc;
typedef long long ll;
typedef long double ld;

const ull millerrabin_prime[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
const ull prime[] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
atomic<ull> cnt(0);
atomic<ull> it(0);
const ull th_count = std::thread::hardware_concurrency();
const ull delta = th_count * 6 - 2;

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

void primelist(vector<ull> *primes, ull n, ull count)
{
    const ull S = 524288;
    ull result = 0, i, j, k, l, start, p;
    Bitslist64 is_prime(S >> 1, true);
    for (k = (++it); k * S <= n; k = (++it)) {
        start = k * S;
        is_prime.reset(true);
        for (l = 0; l < count; l++) {
            p = (*primes)[l];
            for (j = max(((start + p - 1) / p) | 1, p) * p - start; j < S; j += p << 1)
                is_prime.setfalse(j >> 1);
        }
        for (i = 1; i < S and start + i <= n; i += 2) {
            if (is_prime.at(i >> 1)) {
                result++;
                // printf("%u ", i);
            }
        }
    }
    cnt.store((ull)cnt + result);
}
void mbe_sieve(ull n)
{ // multithreading blocked Eratosthenes sieve
    if (n == 1) {
        printf("\nThere are 0 prime numbers in total, and the calculation takes 0.0000 seconds\n");
        return;
    }
    auto start_time = chrono::high_resolution_clock::now(); // 开始计时

    const ull S = 524288;
    vector<ull> primes;
    ull nsqrt = (ull)sqrt(n), count = 0, result = 0, i, j, p;
    Bitslist64 is_prime0((nsqrt >> 1) + 1, true);
    for (i = 3; i <= nsqrt; i += 2) {
        if (is_prime0.at(i >> 1)) {
            primes.push_back(i);
            count++;
            for (j = i * i; j <= nsqrt; j += i << 1)
                is_prime0.setfalse(j >> 1);
        }
    }
    Bitslist64 is_prime(S >> 1, true);
    for (i = 0; i < count; i++) {
        p = primes[i];
        for (j = p * p; j < S; j += p << 1)
            is_prime.setfalse(j >> 1);
    }
    for (i = 3; i < S and i <= n; i += 2) {
        if (is_prime.at(i >> 1)) {
            result++;
            // printf("%u ", i);
        }
    }
    thread *threads = new thread[th_count];
    for (i = 0; i < th_count; i++)
        threads[i] = thread(primelist, &primes, n, count);

    for (i = 0; i < th_count; i++)
        threads[i].join();

    result += (ull)cnt;
    result++; // 把2加进去
    used_time = duration_cast<milliseconds>(high_resolution_clock::now() - start_time).count() / 1000.0;
    printf("\nThere are %llu prime numbers in total, and the calculation takes %.4f seconds\n", result, used_time);
    // delete[] primes;
    delete[] threads;
    cnt = 0;
    it = 0;
}
void PrimeList64(ull start, ull end)
{
    if (start > end) {
        printf("The START cannot be greater than the END!!!\n");
        return;
    }
    if (start == end and start == 1) {
        printf("Within the range of [1, 1], There are 0 prime numbers in total\n");
        return;
    }
    auto start_time = chrono::high_resolution_clock::now(); // 开始计时

    ull i, j, root = (ull)sqrt(end) + 1, count = 1;
    Bitslist64 is_prime((end - 1) / 2, true);
    printf("Within the range of [%llu, %llu], ", start, end);
    start = max(3, start); // 起始值
    if (print_res) {
        printf("there are prime numbers:\n2 ");
        for (i = 3; i < start; i += 2) {
            if (is_prime.at(i >> 1)) {
                for (j = i * i; j <= end; j += i << 1)
                    is_prime.setfalse(j >> 1);
            }
        }
        for (i = start; i < root; i += 2) {
            if (is_prime.at(i >> 1)) {
                count++;
                printf("%llu ", i);
                for (j = i * i; j <= end; j += i << 1)
                    is_prime.setfalse(j >> 1);
            }
        }
        for (i = max(root, start) | 1; i <= end; i += 2) {
            if (is_prime.at(i >> 1)) {
                count++;
                printf("%llu ", i);
            }
        }
    } else {
        for (i = 3; i < start; i += 2) {
            if (is_prime.at(i >> 1)) {
                for (j = i * i; j <= end; j += i << 1)
                    is_prime.setfalse(j >> 1);
            }
        }
        for (i = start; i < root; i += 2) {
            if (is_prime.at(i >> 1)) {
                count++;
                for (j = i * i; j <= end; j += i << 1)
                    is_prime.setfalse(j >> 1);
            }
        }
        for (i = max(root, start) | 1; i <= end; i += 2) {
            if (is_prime.at(i >> 1))
                count++;
        }
    }
    used_time = duration_cast<milliseconds>(high_resolution_clock::now() - start_time).count() / 1000.0;
    printf("\nThere are %llu prime numbers in total, and the calculation takes %.4f seconds\n", count, used_time);
}

bool isprime(ull num)
{
    ull i = 7;
    const ull root = (ull)sqrt(num);
    bool is_p = true;
    if (num == 1) {
        printf("1 is neither a prime nor a composite number\n");
        return false;
    }
    // start_time = chrono::high_resolution_clock::now();

    if (num != 2 and num % 2 == 0)
        is_p = false;
    else if (num != 3 and num % 3 == 0)
        is_p = false;
    else if (num != 5 and num % 5 == 0)
        is_p = false;
    else {
        while (i <= root) {
            if (num % i == 0) {
                is_p = false;
                break;
            }
            i += 4; // 7+4=11
            if (num % i == 0) {
                is_p = false;
                break;
            }
            i += 2; // 11+2=13
            if (num % i == 0) {
                is_p = false;
                break;
            }
            i += 4; // 13+4=17
            if (num % i == 0) {
                is_p = false;
                break;
            }
            i += 2; // 17+2=19
            if (num % i == 0) {
                is_p = false;
                break;
            }
            i += 4; // 19+4=23
            if (num % i == 0) {
                is_p = false;
                break;
            }
            i += 6; // 23+6=29
            if (num % i == 0) {
                is_p = false;
                break;
            }
            i += 2; // 29+2=31
            if (num % i == 0) {
                is_p = false;
                break;
            }
            i += 6; // 31+6=37
        }
    }
    if (is_p) {
        // printf("YES\n");
        return true;
    } else {
        // printf("NO\n");
        return false;
    }
    // used_time = duration_cast<microseconds>(high_resolution_clock::now() - start_time).count() / 1000.0;
    // printf("The calculation takes %.4f milliseconds\n", used_time);
}
inline ull qmul(ull a, ull b, ull mod) { return (__int128)a * b % mod; }

inline ull qpow(ull a, ull n, ull mod)
{
    ull res = 1;
    while (n) {
        if (n & 1)
            res = qmul(res, a, mod);
        a = qmul(a, a, mod);
        n >>= 1;
    }
    return res;
}

inline void output(int mode = 0)
{
    if (mode == 0) {
        used_time = duration_cast<milliseconds>(high_resolution_clock::now() - start_time).count() / 1000.0;
        printf("The calculation takes %.4f seconds\n", used_time);
        return;
    }
    if (mode == 1) {
        used_time = duration_cast<microseconds>(high_resolution_clock::now() - start_time).count() / 1000.0;
        printf("The calculation takes %.4f milliseconds\n", used_time);
    }
}
void Miller_Rabin(ull n) // 判断素数
{
    start_time = chrono::high_resolution_clock::now();
    if (n < 3 or n % 2 == 0) {
        if (n == 2) {
            printf("YES\n");
            output(1);
            return;
        }
        printf("NO\n");
        output(1);
        return;
    }
    ull u = n - 1, t = 0;
    while (u % 2 == 0) {
        u = u >> 1;
        t++;
    }

    for (ull a : millerrabin_prime) {
        ull v = qpow(a, u, n);
        if (v == 1 or v == n - 1 or v == 0)
            continue;
        for (ull j = 1; j <= t; j++) {
            v = qmul(v, v, n);
            if (v == n - 1 and j != t) {
                v = 1;
                break;
            } // 出现一个n-1，后面都是1，直接跳出
            if (v == 1) {
                printf("NO\n");
                output(1);
                return;
            }
        }
        if (v != 1) {
            printf("NO\n");
            output(1);
            return;
        }
    }
    printf("YES\n");
    output(1);
    return;
}

bool check_prime_gmp(const char *n)
{
    mpz_t num, u, base, v, temp;
    mpz_init_set_str(num, n, 10); // num = n
    mpz_init(u);
    mpz_init(base);
    mpz_init(v);
    mpz_init(temp);
    for (ull p : prime) {
        mpz_mod_ui(temp, num, p);
        if (mpz_cmp_ui(temp, 0) == 0) // if(n % p == 0)
            return mpz_cmp_ui(num, p) == 0;
    }
    mpz_sub_ui(u, num, 1); // u = num - 1
    ull t = 0;
    while (mpz_even_p(u)) {  // u % 2 == 0
        mpz_div_ui(u, u, 2); // u = u / 2
        t++;
    }

    // ull ud[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    for (ull a : millerrabin_prime) {
    FLAG:
        mpz_set_ui(base, a);
        mpz_powm(v, base, u, num);                                 // v = qpow(a, u, n);
        if (mpz_cmp_ui(v, 1) == 0 or mpz_cmp(v, u) == 0 or v == 0) // if (v == 1 or v == u or v == 0)
            continue;
        for (ull j = 1; j <= t; j++) {
            mpz_mul(temp, v, v); // v = qmul(v, v, n)
            mpz_mod(v, temp, num);
            if (mpz_cmp(v, u) == 0 and j != t) { // if (v == u and j != t)
                goto FLAG;
                break;
            }
            if (mpz_cmp_ui(v, 1) == 0) // if (v == 1)
                return false;
        }
        if (mpz_cmp_ui(v, 1) != 0) // if (v != 1)
            return false;
    }
    return true;
}
inline bool check_prime(ull n)
{
    if (n == 1) {
        return false;
    }
    for (ull p : prime) {
        if (n % p == 0)
            return n == p;
    }
    ull u = n - 1, t = 0;
    while (u % 2 == 0)
        u /= 2, ++t;
    // ull ud[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    for (ull a : millerrabin_prime) {
        ull v = qpow(a, u, n);
        if (v == 1 or v == n - 1 or v == 0)
            continue;
        for (ull j = 1; j <= t; j++) {
            v = qmul(v, v, n);
            if (v == n - 1 and j != t) {
                v = 1;
                break;
            } // 出现一个n-1，后面都是1，直接跳出
            if (v == 1)
                return 0; // 这里代表前面没有出现n-1这个解，二次检验失败
        }
        if (v != 1)
            return 0; // Fermat检验
    }
    return 1;
}

void a0(ull start, ull end, ull c)
{
    ull count = 0, num;
    for (num = 6 * (start / 6 + c) + 5; num <= end;) {
        if (check_prime(num)) {
            count++;
            printf("%llu ", num);
        }
        num += 2;

        if (num > end)
            break;

        if (check_prime(num)) {
            count++;
            printf("%llu ", num);
        }
        num += delta;
    }
    cnt.store((ull)cnt + count);
}
void a1(ull start, ull end, ull c)
{
    ull count = 0, num;
    for (num = 6 * (start / 6 + c) + 5; num <= end;) {
        if (check_prime(num))
            count++;
        num += 2;

        if (num > end)
            break;

        if (check_prime(num))
            count++;
        num += delta;
    }
    cnt.store((ull)cnt + count);
}
void primelist_multithread(ull start, ull end)
{
    if (start > end) {
        printf("The START cannot be greater than the END!!!\n");
        return;
    }
    printf("Within the range of [%llu, %llu], there are prime numbers:\n", start, end);
    start_time = chrono::high_resolution_clock::now();
    if (start <= 2 and 2 <= end) {
        cnt++;
    }
    if (start <= 3) {
        if (3 <= end)
            cnt++;
        else {
            printf("\nThere are %llu prime numbers in total", (ull)cnt);
            output();
            return;
        }
    }
    thread *threads = new thread[th_count];
    ull i;
    if (print_res) {
        if (start % 6 == 0 or start % 6 == 1) {
            if (start != 1 and check_prime(start)) {
                printf("%llu ", start);
                cnt++;
            }
        }
        for (i = 0; i < th_count; i++)
            threads[i] = thread(a0, start, end, i);
    } else {
        if (start % 6 == 0 or start % 6 == 1) {
            if (start != 1 and check_prime(start))
                cnt++;
        }
        for (i = 0; i < th_count; i++)
            threads[i] = thread(a1, start, end, i);
    }
    for (i = 0; i < th_count; i++)
        threads[i].join();

    used_time = duration_cast<milliseconds>(high_resolution_clock::now() - start_time).count() / 1000.0;
    printf("\nThere are %llu prime numbers in total, and the calculation takes %.4f seconds\n", (ull)cnt, used_time);
    cnt = 0;
    delete[] threads;
}

void factorize(ull n)
{
    if (n == 1) {
        printf("Invalid input!!!\n");
        return;
    }
    ull root = (ull)sqrt(n), count = 0, types = 0, indices = 1, i = 7;
    vector<ull> factors;
    printf("Factors         Exponent\n");
    start_time = chrono::high_resolution_clock::now();
    while (n != 1 && 2 <= root) {
        if (n % 2 == 0) {
            n = n / 2;
            factors.push_back(2);
            count++;
        } else
            break;
    }
    while (n != 1 && 3 <= root) {
        if (n % 3 == 0) {
            n = n / 3;
            factors.push_back(3);
            count++;
        } else
            break;
    }
    while (n != 1 && 5 <= root) {
        if (n % 5 == 0) {
            n = n / 5;
            factors.push_back(5);
            count++;
        } else
            break;
    }
    while (n != 1 && i <= root) {
        if (n % i == 0) {
            n = n / i;
            factors.push_back(i);
            count++;
        } else
            i += 4; // 7+4=11
        if (n % i == 0) {
            n = n / i;
            factors.push_back(i);
            count++;
        } else
            i += 2; // 11+2=13
        if (n % i == 0) {
            n = n / i;
            factors.push_back(i);
            count++;
        } else
            i += 4; // 13+4=17
        if (n % i == 0) {
            n = n / i;
            factors.push_back(i);
            count++;
        } else
            i += 2; // 17+2=19
        if (n % i == 0) {
            n = n / i;
            factors.push_back(i);
            count++;
        } else
            i += 4; // 19+4=23
        if (n % i == 0) {
            n = n / i;
            factors.push_back(i);
            count++;
        } else
            i += 6; // 23+6=29
        if (n % i == 0) {
            n = n / i;
            factors.push_back(i);
            count++;
        } else
            i += 2; // 29+2=31
        if (n % i == 0) {
            n = n / i;
            factors.push_back(i);
            count++;
        } else
            i += 6; // 31+6=37
    }
    if (n != 1) {
        factors.push_back(n);
        count++;
    }
    factors.push_back(0);
    for (i = 0; i < count; i++) {
        if (factors[i] != factors[(ull)i + 1]) {
            cout << factors[i] << "		" << indices << endl;
            indices = 1;
            types++;
        } else
            indices++;
    }
    used_time = duration_cast<microseconds>(high_resolution_clock::now() - start_time).count() / 1000.0;
    printf("The calculation takes %.4f milliseconds\n", used_time);
}

bool execute_command(string command)
{
    if (command == "/m") {
        cout << mode_info << choose_mode_info;
        return true;
    }
    if (command == "/h") {
        cout << help_info;
        return false;
    }
    if (command == "/r") {
        print_res = !print_res;
        return false;
    }
    if (command == "/q") {
        exit(0);
    }
    printf("Invalid command!!!\n");
    return false;
}

bool StringtoUll(ull *arg, string str)
{
    ull num = 0;
    if (str == "" or str.size() > 20)
        return false;
    for (uint i = 0; i < str.size(); i++) {
        if (str.at(i) < '0' or str.at(i) > '9')
            return false;
    }
    for (uint i = 0; i < str.size(); i++)
        num = num * 10 + str.at(i) - '0';

    if (num >= 1e19 or num <= 0)
        return false;
    *arg = (ull)num;
    return true;
}

ull get_argument(string message = "")
{
    ull arg = 0;
    string entry;

    while (1) {
        cout << message;
        getline(cin, entry);
        if (entry[0] == '/') {
            restart = execute_command(entry);
            if (restart)
                break; // 重新选择
        } else if (StringtoUll(&arg, entry))
            break; // 转成数字后返回

        else
            printf("Invalid input!!!\n");
    }
    return arg;
}

string get_argument_str(string message = "")
{
    string entry;

    while (1) {
        cout << message;
        getline(cin, entry);
        if (entry[0] == '/') {
            restart = execute_command(entry);
            if (restart)
                break; // 重新选择
        } else {
            for (uint i = 0; i < entry.size(); i++) {
                if (entry.at(i) < '0' or entry.at(i) > '9')
                    printf("Invalid input!!!\n");

                else
                    return entry;
            }
        }
    }
    return "";
}

int main()
{
    string entry;
    printf("============================\n=Welcome to this calculator=\n============================\n");
    cout << help_info << choose_mode_info << mode_info;
    while (1) {
        bool t1;
        restart = false;
        getline(cin, entry);
        if (entry[0] == '/') {
            restart = execute_command(entry);
            if (restart)
                continue;                 // 重新选择
        } else if (entry.length() == 1) { // 如果字符串长度为1
            ull arg1 = 0, arg2 = 0;
            UPmath::BigInteger heihiehie;
            switch (entry[0]) {
            case '1': // 模式1，生成质数序列
                printf("[Current mode]: Generate prime number list-1\n");
                arg1 = get_argument("Please enter the starting value: ");
                if (restart)
                    break; // 重新选择

                arg2 = get_argument("Please enter the end value: ");
                if (restart)
                    break; // 重新选择

                PrimeList64(arg1, arg2);
                cout << choose_mode_info;
                break;

            case '2': // 模式2，判断一个数是否为质数
                printf("[Current mode]: Primality test\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择
                start_time = high_resolution_clock::now();
                t1 = isprime(arg1);
                output(1);
                if (t1)
                    printf("YES\n");

                else
                    printf("YES\n");

                cout << choose_mode_info;
                break;

            case '3': // 模式3，分解质因数
                printf("[Current mode]: Factorization\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break;

                factorize(arg1);
                cout << choose_mode_info;
                break;

            case '4': // 模式4，生成质数序列2
                printf("[Current mode]: Generate prime number list-2\n");
                arg1 = get_argument("Please enter the starting value: ");
                if (restart)
                    break; // 重新选择

                arg2 = get_argument("Please enter the end value: ");
                if (restart)
                    break; // 重新选择

                primelist_multithread(arg1, arg2);
                cout << choose_mode_info;
                break;

            case '5': // 模式5 分块埃筛
                printf("[Current mode]: Multi-threading blocked Eratosthenes sieve\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择
                mbe_sieve(arg1);
                cout << choose_mode_info;
                break;

            case '6': // 模式6，米勒-卡宾素性检测
                printf("[Current mode]: Miller-Rabin \n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择

                Miller_Rabin(arg1);
                cout << choose_mode_info;
                break;

            case '7': // mode-7, debug-1
                printf("[Current mode]: Debug-1\n");
                entry = get_argument_str("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择
                start_time = high_resolution_clock::now();
                t1 = check_prime_gmp(entry.c_str());
                output(1);
                if (t1)
                    printf("YES\n");
                else
                    printf("NO\n");
                cout << choose_mode_info;
                break;

            case '8': // mode-8，debug-2
                printf("[Current mode]: Debug-2\n");
                entry = get_argument_str("Please enter a positive integer: ");

                if (restart)
                    break; // 重新选择
                start_time = high_resolution_clock::now();
                heihiehie = UPmath::BigInteger(entry);
                t1 = heihiehie.weakerBailliePSWPrimeTest();
                output(1);
                if (t1)
                    printf("YES\n");
                else
                    printf("NO\n");

                cout << choose_mode_info;
                break;

            default:
                cout << invaild_mode;
                break;
            }
        } else
            cout << invaild_mode;
    }
    return 0;
}