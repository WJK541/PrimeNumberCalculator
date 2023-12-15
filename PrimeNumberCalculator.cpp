#include "head.h"

const string invaild_mode = "Invalid mode!!!\n";
const string help_info = "\n[Commands]:\n/h  Get commands\n/m  Re-select mode\n/q  Quit the program\n";
const string mode_info = "1)Generate prime number list-1    2)Primality test    3)Factorization    4)Generate prime number list-2\n5)Eratosthenes sieve    6)Miller-Rabin    7)Debug-1    8)Debug-2\n";
const string choose_mode_info = "\n[Select mode]:";

bool execute_command(string cmd)
{
    if (commands.find(cmd) != commands.end()) {
        bool (*funcPointer)() = commands[cmd];
        return funcPointer();
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

int main()
{
    commands["/m"] = []() -> bool {  cout << mode_info << choose_mode_info;return true; };
    commands["/h"] = []() -> bool {  cout << help_info;return false; };
    commands["/r"] = []() -> bool {  print_res = !print_res;return false; };
    commands["/q"] = []() -> bool {  exit(0);return true; };
    cout << "============================\n"
         << "=Welcome to this calculator=\n"
         << "============================\n";
    cout << help_info << choose_mode_info << mode_info;
    while (1) {
        string entry;
        bool t1, t2;
        getline(cin, entry);
        if (entry[0] == '/') {
            restart = execute_command(entry);
            if (restart)
                continue;
        } else if (entry.length() == 1) {
            ull arg1 = 0, arg2 = 0;
            vector<ull> v1 = {};
            vector<vector<ull> *> v2;
            switch (entry[0]) {
            case '1': // 模式1，生成质数序列
                printf("[Current mode]: Generate prime number list-1\n");
                arg1 = get_argument("Please enter the starting value: ");
                if (restart)
                    break; // 重新选择

                arg2 = get_argument("Please enter the end value: ");
                if (restart)
                    break; // 重新选择

                START_TIMING
                v1 = PrimeList64(arg1, arg2);
                END_TIMING_MS

                if (print_res) {
                    for (int i = 0; i < v1.size(); i++) {
                        cout << v1[i] << ' ';
                    }
                }
                cout << '\n'
                     << v1.size() << "\nFinish\n";
                cout << choose_mode_info;
                break;

            case '2': // 模式2，判断一个数是否为质数
                      // 计划删除
                printf("[Current mode]: Primality test\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择
                t1 = MillerRabin(arg1);
                output(t1);

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
                arg1 = get_argument("Please enter the start value: ");
                if (restart)
                    break; // 重新选择

                arg2 = get_argument("Please enter the end value: ");
                if (restart)
                    break; // 重新选择

                START_TIMING
                v2 = PrimeListMT(arg1, arg2);
                END_TIMING_MS
                output(v2);
                
                cout << choose_mode_info;
                break;

            case '5': // 模式5 分块埃筛
                printf("[Current mode]: Multi-threaded blocked Eratosthenes sieve\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择
                START_TIMING
                v2 = EratosthenesSieve(arg1);
                END_TIMING_MS
                output(v2);

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
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择
                START_TIMING
                // v1 = EratosthenesSieve(arg1);
                END_TIMING_MS

                if (print_res) {
                    for (ull i : v1) {
                        cout << i << ' ';
                    }
                }

                cout << '\n'
                     << v1.size() << "\nFinish\n";
                cout << choose_mode_info;
                break;

            case '8': // mode-8，debug-2
                printf("[Current mode]: Debug-2\n");
                arg1 = get_argument("Please enter the start value: ");
                if (restart)
                    break; // 重新选择

                arg2 = get_argument("Please enter the end value: ");
                if (restart)
                    break; // 重新选择
                cout << '\n'
                     << 1 << '\n';
                START_TIMING
                for (int i = 0; i < 10; i++) {
                    PrimeListMT(arg1, arg2);
                }
                END_TIMING_MS

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