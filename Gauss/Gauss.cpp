#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const double EPS = 0.000001;
const double INF = 1e9 + 7;

int gauss(vector<vector<double>> A, vector<double>& x)
{
    int n = (int)A.size();
    int m = (int)A[0].size()-1;

    vector<int> where(m, -1);
    for (int col = 0, row = 0; col < m && row < n; ++col)
    {
        int sel = row;
        for (int i = row; i < n; ++i)
        {
            if (abs(A[i][col]) > abs(A[sel][col]))
                sel = i;
        }

        if (abs(A[sel][col]) < EPS)
            continue;

        for (int i = col; i <= m; ++i)
        {
            swap(A[sel][i], A[row][i]);
        }
        where[col] = row;

        for (int i = 0; i < n; ++i)
        {
            if (i != row)
            {
                double c = A[i][col] / A[row][col];
                for (int j = col; j <= m; ++j)
                {
                    A[i][j] -= A[row][j] * c;
                }
            }
        }

        ++row;
    }

    for (int i = 0; i < m; ++i)
    {
        if (where[i] != -1)
            x[i] = A[where[i]][m] / A[where[i]][i];
    }
        
    for (int i = 0; i < n; ++i) 
    {
        double sum = 0;
        for (int j = 0; j < m; ++j) 
        {
            sum += x[j] * A[i][j];
        }
            
        if (abs(sum - A[i][m]) > EPS)
            return 0;
    }

    for (int i = 0; i < m; ++i) 
    {
        if (where[i] == -1)
            return INF;
    }

    return 1;
}

int main()
{
    setlocale(LC_ALL, "Russian");

    // Ввод данных
    int n, m;
    cout << "Введите количество уравнений: ";
    cin >> n;
    cout << endl;
    cout << "Введите количество неизвестных: ";
    cin >> m;
    cout << endl;

    vector<vector<double>> A(n, vector<double>(m+1));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << "A[" << i << "][" << j << "] = ";
            cin >> A[i][j];
        }
    }
    cout << endl;

    for (int i = 0; i < n; i++)
    {
        cout << "b[" << i << "] = ";
        cin >> A[i][m];
    }
    cout << endl;

    vector<double> x(m, 0);
    
    int ans = gauss(A, x);

    // Вывод ответа

    if (ans == 1) 
    {
        cout << "Решение системы:" << endl;
        for (int i = 0; i < m; i++)
        {
            cout << "x[" << i << "] = " << x[i] << endl;
        }
    }
    else if (ans == INF) 
    {
        cout << "Система имеет бесконечно много решений" << endl;
    }
    else if(ans == 0)
    {
        cout << "Система не имеет решений" << endl;
    }

    return 0;
}