#include <vector>
#include <iostream>
#include <omp.h>
#include <windows.h>

using namespace std;

int64_t inverse(int64_t a, int64_t b)
{
    int64_t b0 = b, t, q;
    int64_t x0 = 0, x1 = 1;
    if (b == 1) return 1;
    while (a > 1) {
        q = a / b;
        t = b, b = a % b, a = t;
        t = x0, x0 = x1 - q * x0, x1 = t;
    }
    if (x1 < 0) x1 += b0;
    return x1;
}

const int64_t mod = 180143985094819841LL;
const int64_t root = 7776LL;
const int64_t root_1 = inverse(root, mod);
const int64_t two_1 = inverse(2, mod);
const int64_t root_pw = 1LL << 55;

int log2(int64_t n)
{
    int64_t p; int r;
    for (p = 1, r = 0; p < n; p <<= 1, r++);
    return r;
}

int log10(int n)
{
    if (n / 10000000000) return 10;
    if (n / 1000000000) return 9;
    if (n / 100000000) return 8;
    if (n / 10000000) return 7;
    if (n / 1000000) return 6;
    if (n / 100000) return 5;
    if (n / 10000) return 4;
    if (n / 1000) return 3;
    if (n / 100) return 2;
    if (n / 10) return 1;
    return 0;
}

int64_t pow10(int n)
{
    if (n == 10) return 10000000000;
    if (n == 9) return 1000000000;
    if (n == 8) return 100000000;
    if (n == 7) return 10000000;
    if (n == 6) return 1000000;
    if (n == 5) return 100000;
    if (n == 4) return 10000;
    if (n == 3) return 1000;
    if (n == 2) return 100;
    if (n == 1) return 10;
    if (n == 0) return 1;
    return 0;
}

inline int64_t add_mod(int64_t a, int64_t b)
{
    int64_t r = a + b;
    r -= mod;
    r += (r >> 63) & mod;
    return r;
}

inline int64_t sub_mod(int64_t a, int64_t b)
{
    int64_t r = a - b;
    r += (r >> 63) & mod;
    return r;
}

inline int64_t mod_mul(int64_t a, int64_t b)
{
    int64_t r, rh, t;
    asm (" \
         imulq %3; \
         movq %%rax, %%rcx; \
         shrq $55, %%rax; \
         shlq $9, %%rdx; \
         orq %%rdx, %%rax; \
         mulq %4; \
         shrq $2, %%rdx; "
         : "=d"(rh), "=c"(t)
         : "%a"(a), "r"(b), "r"(0xCCCCCCCCCCCCCCCD), "r"(mod));
    rh *= mod;
    r = t - rh;
    return r < 0 ? r + mod : r;
}

int64_t mod_mul_two_1(int64_t a)
{
    int64_t r = two_1 + (a >> 1);
    int64_t d = (a & 1 == 1 ? mod : two_1);
    return r >= d ? r - d : r;
}

inline void ntt_2_point(vector<int64_t> &a, int64_t i, bool invert)
{
    int64_t u = a[i];
    int64_t v = a[i + 1];
    a[i] = add_mod(u, v);
    a[i + 1] = sub_mod(u, v);
    if (invert) {
        a[i] = mod_mul_two_1(a[i]);
        a[i + 1] = mod_mul_two_1(a[i + 1]);
    }
}

inline void ntt_4_point(vector<int64_t> &a, int64_t i, int64_t wlen, bool invert)
{
    int64_t u, v, t;
    u = add_mod(a[i], a[i + 2]);
    v = sub_mod(a[i], a[i + 2]);
    a[i] = u;
    a[i + 2] = v;

    t = mod_mul(a[i + 3], wlen);
    u = add_mod(a[i + 1], t);
    v = sub_mod(a[i + 1], t);
    a[i + 1] = u;
    a[i + 3] = v;

    if (invert) {
        a[i] = mod_mul_two_1(a[i]);
        a[i + 1] = mod_mul_two_1(a[i + 1]);
        a[i + 2] = mod_mul_two_1(a[i + 2]);
        a[i + 3] = mod_mul_two_1(a[i + 3]);
    }
}

void ntt(vector<int64_t> &a, bool invert, int p)
{
    int64_t n = a.size();

    #pragma omp parallel for if(p) schedule(static)
    for (int64_t i = 1; i < n; i++) {
        uint64_t j = i;
        j = (((j & 0xaaaaaaaa) >> 1) | ((j & 0x55555555) << 1));
        j = (((j & 0xcccccccc) >> 2) | ((j & 0x33333333) << 2));
        j = (((j & 0xf0f0f0f0) >> 4) | ((j & 0x0f0f0f0f) << 4));
        j = (((j & 0xff00ff00) >> 8) | ((j & 0x00ff00ff) << 8));
        j = (((j & 0xffff0000) >> 16) | ((j & 0x0000ffff) << 16));
        j = ((j >> 32) | (j << 32)) >> (64 - log2(n));
        if (j > i)
            swap(a[i], a[j]);
    }

    for (int64_t len = 2; len <= n; len <<= 1) {
        int64_t wlen = invert ? root_1 : root;
        for (int64_t i = len; i < root_pw; i <<= 1)
            wlen = mod_mul(wlen, wlen);
        #pragma omp parallel for if(p) schedule(static)
        for (int64_t i = 0; i < n; i += len) {
            int64_t w = 1;
            if (len == 2)
                ntt_2_point(a, i, invert);
            else if (len == 4)
                ntt_4_point(a, i, wlen, invert);
            else
                for (int64_t j = 0; j < len / 2; j++) {
                    int64_t u = a[i + j], v = mod_mul(a[i + j + len / 2], w);
                    a[i + j] = add_mod(u, v);
                    a[i + j + len / 2] = sub_mod(u, v);
                    if (invert) {
                        a[i + j] = mod_mul_two_1(a[i + j]);
                        a[i + j + len / 2] = mod_mul_two_1(a[i + j + len / 2]);
                    }
                    w = mod_mul(w, wlen);
                }
        }
    }
}

void ntt_mul(vector<int> &a, vector<int> &b, int p)
{
    int64_t n = 1;
    while (n < a.size() + b.size())
        n <<= 1;

    vector<int64_t> x(a.begin(), a.end());
    x.resize(n, 0);
    vector<int64_t> y(b.begin(), b.end());
    y.resize(n, 0);

    ntt(x, false, p);
    ntt(y, false, p);
    #pragma omp parallel for if(p)
    for (int64_t i = 0; i < n; i++)
        x[i] = mod_mul(x[i], y[i]);
    ntt(x, true, p);

    a.resize(n, 0);
    int64_t r = 0, w;
    for (int64_t i = 0; i < n; i++) {
        w = x[i] + r;
        a[i] = w % 10000;
        r = w / 10000;
    }
}

void print(vector<int> a)
{
    auto e = a.end() - 1;
    while (!*e) e--;

    for (int i = 4; i > 0; i--)
        if(*e % pow10(i) / pow10(i - 1))
        {
            for (int j = i; j > 0; j--)
                cout << *e % pow10(j) / pow10(j - 1);
            break;
        }

    for (auto i = --e; i >= a.begin(); i--)
        for (int j = 4; j > 0; j--)
            cout << *i % pow10(j) / pow10(j - 1);
}

void factorial(int n)
{
    int64_t i, j, r, w, t;

    int thread_count = 28; // omp_get_max_threads();
    while (thread_count & thread_count - 1)
        thread_count = thread_count & thread_count - 1;

    int parts_count = thread_count;
    while (parts_count < n / 10000)
        parts_count <<= 1;
    int k = parts_count / thread_count;
    int m = 0;

    vector<vector<int>> a(parts_count);
    vector<int> l(parts_count, 1);

    for (i = 0; i < parts_count; i++)
    {
        a[i].resize(n * log10(n) / parts_count, 0);
        a[i][0] = 1;
    }

    omp_set_num_threads(thread_count);
    #pragma omp parallel for private(t, i, j, r, w) firstprivate(m) schedule(dynamic)
    for (j = 1; j <= n; j++)
    {
        t = omp_get_thread_num() * k + m;
        i = 0;
        r = 0;
        while (i < l[t] || r > 0)
        {
            w = (int64_t)a[t][i] * j + r;
            r = w / 100000000;
            a[t][i] = w % 100000000;
            if (i == l[t]) l[t]++;
            i++;
        }
        if(++m == k) m = 0;
    }

    #pragma omp parallel for private(i, j)
    for (i = 0; i < parts_count; i++)
    {
        vector<int> b;
        j = 0;
        while (j < l[i])
        {
            b.push_back(a[i][j] % 10000);
            b.push_back(a[i][j] / 10000);
            j++;
        }
        a[i] = b;
    }

    for (i = parts_count / 2; i > thread_count / 8; i /= 2)
    {
        #pragma omp parallel for
        for (j = 0; j < i; j++)
            ntt_mul(a[j], a[j + i], 0);
    }
    for (i = thread_count / 8; i > 0; i /= 2)
    {
        for (j = 0; j < i; j++)
            ntt_mul(a[j], a[j + i], 1);
    }
    print(a[0]);
}

int main(int argc, char** argv)
{
    DWORD t = GetTickCount();
    factorial(atoi(argv[1]));
    cerr << GetTickCount() - t << endl;
}
