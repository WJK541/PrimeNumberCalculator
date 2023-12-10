#include "head.h"

const string invaild_mode = "Invalid mode!!!\n";
const string help_info = "\n[Commands]:\n/h  Get commands\n/m  Re-select mode\n/r  Turn on/off output result\n/q  Quit the program\n";
const string mode_info = "1)Generate prime number list-1    2)Primality test    3)Factorization    4)Generate prime number list-2\n5)Eratosthenes sieve    6)Miller-Rabin    7)Debug-1    8)Debug-2\n";
const string choose_mode_info = "\n[Select mode]:";
system_clock::time_point start_time;
double elapsed_time;
bool restart = false;
bool print_res = false;

atomic<ull> it(0);
atomic<ull> cnt(0);
mpz_t ZERO, ONE, TWO;

void primelist(vector<ull> *primes, ull n, ull count)
{
    const ull S = 524288;
    ull i, j, k, l, start, p;
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
                primes->push_back(i);
            }
        }
    }
}
vector<ull> MTBEratosthenesSieve(ull n)
{ // multithreaded blocked Eratosthenes sieve
    const ull S = 524288;
    vector<ull> primes;
    ull nsqrt = (ull)sqrt(n), count = 0, i, j, p;
    Bitslist64 is_prime0((nsqrt >> 1) + 1, true);
    for (i = 3; i <= nsqrt; i += 2) {
        if (is_prime0.at(i >> 1)) {
            primes.push_back(i);
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
            primes.push_back(i);
            // printf("%u ", i);
        }
    }
    vector<thread> threads;
    for (i = 0; i < th_count; i++)
        threads.push_back(thread(primelist, &primes, n, count));

    for (i = 0; i < th_count; i++)
        threads[i].join();

    cnt = 0;
    it = 0;
    return primes;
}

inline bool millerrabin(const mpz_t n, const mpz_t nm1, mpz_t exp, mpz_t p, uint t)
{
    mpz_powm(p, p, exp, n); // p = (p ^ exp) % n

    if (mpz_cmp(p, ONE) == 0 or mpz_cmp(p, nm1) == 0) {
        return true;
    }
    mpz_t temp;
    mpz_init(temp);
    while (--t > 0) {
        mpz_mul(temp, p, p); // p = (p * p) % n
        mpz_mod(p, temp, n);
        if (mpz_cmp(p, nm1) == 0) { // p == n - 1
            mpz_clear(temp);
            return true;
        }
    }
    mpz_clear(temp);
    return false;
}

bool isprime_gmp(ll num)
{
    bool isp = true;
    mpz_t n, p, temp;
    mpz_init_set_ui(n, num); // num = n

    mpz_init(p);
    mpz_init(temp);

    if (mpz_cmp_ui(n, 3) < 0 or mpz_even_p(n) == 1) {
        return (mpz_cmp(n, TWO) == 0);
    }

    for (ull p1 : prime) {
        mpz_set_ui(p, p1);
        mpz_mod(temp, n, p);
        if (mpz_cmp(temp, ZERO) == 0) {
            mpz_clear(temp);
            return (mpz_cmp(n, p) == 0);
        }
    }

    mpz_t nm1, m;
    mpz_init(nm1);
    mpz_init(m);
    mpz_sub(nm1, n, ONE); // nm1 = n - 1
    mpz_set(m, nm1);
    uint e = 0;
    while (mpz_even_p(m)) {    // u % 2 == 0
        mpz_cdiv_q(m, m, TWO); // u = u / 2
        e++;
    }
    for (ull p2 : millerrabin_prime) {
        mpz_set_ui(p, p2);
        isp = isp && millerrabin(n, nm1, m, p, e);
    }
    return isp;
}

int _JacobiSymbolImpl(ll numerator, ll denominator);
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
    if (__gcd(numerator, denominator) != 1)
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
        g = __gcd(ds, n);
        if (g != 1 or g != -1)
            return 0;

        if (JacobiSymbol(ds, n) == -1)
            return ds;

        D += 2;
        s = -s;
    }
}
ll getD_debug(ll n)
{
    ll D = 5;
    ll s = 1;
    ll ds, g;
    int j;
    for (;;) {
        ds = D * s;
        g = __gcd(ds, n);
        if (g != 1 or g != -1)
            return 0;

        j = JacobiSymbol(ds, n);
        cout << D << ": " << j << '\n';
        if (j == -1)
            return ds;

        D += 2;
        s = -s;
    }
}
bool lucas(ll n)
{ /*
     for (ll p : prime) {
         if (n % p == 0)
             return n == p;
     }
    ll D = 5;
    ll s = 1;
    ll ds, g;
    for (;;) {
        ds = D * s;
        g = __gcd(ds, n);
        if (abs(g) != 1)
            return 0;

        if (JacobiSymbol(ds, n) == -1)
            break;

        D += 2;
        s = -s;
    }*/
    return lucassequence(n, getD_debug(n));
}

/*
Return 2: N must be a prime number
Return 1: N probably be a prime number
Return 0: N must not be a prime number
*/
int strongPrimalityTest(ull n)
{
    for (ull p : prime) {
        if (n % p == 0)
            if (n == p)
                return 2;
            else
                return 0;
    }
    // n - 1 = r * 2 ^ s
    ull r = (n - 1) >> 1;
    uint s = 1;
    while (r % 2 == 0) {

        r >>= 1;
        s++;
    }
    ull base = 2;
    ull p = 1;
    while (r) {
        if (r & 1)
            p = qmul(p, base, n);
        base = qmul(base, base, n);
        r >>= 1;
    }
    // p = (p ^ r) % n

    ull nMinus1 = n - 1;
    if (p == 1 or p == nMinus1) {
        return lucas(n) ? 1 : 0;
    }
    while (--s > 0) {
        p = qmul(p, p, n);
        // p = (p * p) % n
        if (p == nMinus1) {
            return lucas(n) ? 1 : 0;
        }
    }
    return 0;
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
        } else if (test_set.find(entry) != test_set.end()) {
            return test_set[entry];
        } else {
            for (uint i = 0; i < entry.size(); i++) {
                if (entry.at(i) < '0' or entry.at(i) > '9')
                    printf("Invalid input???\n");
                else
                    return entry;
            }
        }
    }
    return "";
}

bool gmp_miller(ll num)
{
    mpz_t n;
    mpz_init_set_si(n, num);
    return mpz_millerrabin(n, 30) > 0;
}
bool gmp_miller_str(string num)
{
    mpz_t n;
    mpz_init_set_str(n, num.c_str(), 10);
    return mpz_millerrabin(n, 30) > 0;
}
ll gmp_jacobi(ll d, ll n)
{
    mpz_t D, N;
    mpz_init_set_si(D, d);
    mpz_init_set_si(N, n);
    return mpz_jacobi(D, N);
}

int main()
{
    cout << __gcd(-5, 5719) << '\n';
    string entry;
    mpz_init_set_ui(ZERO, 0);
    mpz_init_set_ui(ONE, 1);
    mpz_init_set_ui(TWO, 2);
    printf("============================\n=Welcome to this calculator=\n============================\n");
    cout << help_info << choose_mode_info << mode_info;
    while (1) {
        bool t1, t2;
        restart = false;
        getline(cin, entry);
        if (entry[0] == '/') {
            restart = execute_command(entry);
            if (restart)
                continue;                 // 重新选择
        } else if (entry.length() == 1) { // 如果字符串长度为1
            ull arg1 = 0, arg2 = 0;
            switch (entry[0]) {
            case '1': // 模式1，生成质数序列
                printf("[Current mode]: Generate prime number list-1\n");
                arg1 = get_argument("Please enter the starting value: ");
                if (restart)
                    break; // 重新选择

                arg2 = get_argument("Please enter the end value: ");
                if (restart)
                    break; // 重新选择

                for (ull p : PrimeList64(arg1, arg2)) {
                    cout << p << ' ';
                }
                cout << "\nFinish\n";
                cout << choose_mode_info;
                break;

            case '2': // 模式2，判断一个数是否为质数
                      // 计划删除
                printf("[Current mode]: Primality test\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择
                t1 = MillerRabin(arg1);
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

                for (auto i : PrimeFactorization(arg1)) {
                    cout << (i.first) << "   " << (i.second) << '\n';
                }
                cout << "\nFinish\n";
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

                for (auto i : PrimeListMT(arg1, arg2)) {
                    cout << i << ' ';
                }
                cout << "\nFinish\n";
                cout << choose_mode_info;
                break;

            case '5': // 模式5 分块埃筛
                printf("[Current mode]: Multi-threading blocked Eratosthenes sieve\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择

                for (ull i : MTBEratosthenesSieve(arg1)) {
                    cout << i << ' ';
                }
                cout << "\nFinish\n";
                cout << choose_mode_info;
                break;

            case '6': // 模式6，米勒-卡宾素性检测
                printf("[Current mode]: Miller-Rabin\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择

                t1 = MillerRabin(arg1);
                if (t1)
                    cout << "YES";
                else
                    cout << "NO";

                cout << choose_mode_info;
                break;

            case '7': // mode-7, debug-1
                printf("[Current mode]: Debug-1\n");
                arg1 = get_argument("Please enter the starting value: ");
                if (restart)
                    break; // 重新选择

                arg2 = get_argument("Please enter the end value: ");
                if (restart)
                    break; // 重新选择

                for (arg1 = arg1 | 1; arg1 < arg2; arg1 += 2) {
                    t1 = ((strongPrimalityTest(arg1) > 0) == gmp_miller(arg1));
                    if (not t1)
                        cout << arg1 << '\n';
                }

                cout << choose_mode_info;
                break;

            case '8': // mode-8，debug-2
                printf("[Current mode]: Debug-2\n");
                arg1 = get_argument("Please enter a positive integer: ");

                if (restart)
                    break; // 重新选择
                start_time = high_resolution_clock::now();
                t1 = lucas(arg1);
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