#include "head.h"

int main()
{
    //read();
    
    cout << "============================\n"
         << "=Welcome to this calculator=\n"
         << "============================\n"
         << help_info << mode_info << choose_mode_info;
    while (1) {
        string entry;
        getline(cin, entry);
        if (entry[0] == '/') {
            if (execute_command(entry))
                continue;
        } else if (entry.length() == 1) {
            ull arg1 = 0, arg2 = 0;

            // Result r = Result(false);
            switch (entry[0]) {
            case '1': // 模式1，生成质数序列
                printf("[Current mode]: Generate prime number list-1\n");
                arg1 = get_argument("Please enter the starting value: ");
                if (restart)
                    break; // 重新选择

                arg2 = get_argument("Please enter the end value: ");
                if (restart)
                    break; // 重新选择

                {
                    auto r = Result(PrimeList64(arg1, arg2));
                    r.output();
                    r.write();
                    verify(r.var);
                }

                cout << choose_mode_info;
                break;

            case '2': // 模式2，判断一个数是否为质数
                      // 计划删除
                printf("[Current mode]: Primality test\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择
                {
                    auto r = Result(MillerRabin(arg1));
                    r.output();
                }

                cout << choose_mode_info;
                break;

            case '3': // 模式3，分解质因数
                printf("[Current mode]: Factorization\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break;

                {
                    auto r = Result(PrimeFactorization(arg1));
                    r.output();
                }

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

                {
                    auto r = Result(PrimeListMT(arg1, arg2));
                    r.output();
                    r.write();
                    verify(r.var);
                }

                cout << choose_mode_info;
                break;

            case '5': // 模式5 分块埃筛
                printf("[Current mode]: Multi-threaded blocked Eratosthenes sieve\n");
                arg1 = get_argument("Please enter the starting value: ");
                if (restart)
                    break; // 重新选择

                arg2 = get_argument("Please enter the end value: ");
                if (restart)
                    break; // 重新选择

                {
                    auto r = Result(EratosthenesSieve(arg1, arg2));
                    r.output();
                    //r.write();
                    verify(r.var);
                }

                cout << choose_mode_info;
                break;

            case '6': // 模式6，米勒-卡宾素性检测
                printf("[Current mode]: Miller-Rabin\n");
                arg1 = get_argument("Please enter a positive integer: ");
                if (restart)
                    break; // 重新选择

                {
                    auto r = Result(MillerRabin(arg1));
                    r.output();
                }

                cout << choose_mode_info;
                break;

            case '7': // mode-7, debug-1
                printf("[Current mode]: Debug-1\n");
                cout << "nothing\n";
                cout << choose_mode_info;
                break;

            case '8': // mode-8，debug-2
                printf("[Current mode]: Debug-2\n");
                cout << "nothing\n";
                cout << choose_mode_info;
                break;

            default:
                cout << invalid_mode;
                break;
            }
        } else
            cout << invalid_mode;
    }
    inputFile.close();
    return 0;
}