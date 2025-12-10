#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime> 

using namespace std;

// Параметры задачи //

// Область D: вершины треугольника
const double Cx = -3.0, Cy = 0.0;
const double Ax = 3.0, Ay = 0.0;
const double Bx = 0.0, By = 4.0;

// Вершины прямоугольника П
const double A1 = -4.0, B1 = 4.0; // x-границы прямоугольника
const double A2 = -1.0, B2 = 5.0;  // y-границы прямоугольника

// Параметры сетки
const int M = 400, N = 600;
const double h1 = (B1 - A1) / M;
const double h2 = (B2 - A2) / N;

const double h = max(h1, h2);
const double eps = h*h;    // малый параметр метода фиктивных областей
const double delta = 1e-8;      // точность итерационного метода

// Координата x пересечения с прямых
double get_x_inter(double y, double Ax, double Ay, double Bx, double By) {
    return Bx + ((Ax - Bx) / (Ay - By))  * (y - By);
}

// Координата y пересечения с прямых
double get_y_inter(double x, double Ax, double Ay, double Bx, double By) {
    return Ay + ((By - Ay) / (Bx - Ax)) * (x - Ax);
}

// Проверка, лежит ли точка правее прямой
bool is_right(double x, double y, double Ax, double Ay, double Bx, double By){
    double x_inter = get_x_inter(y, Ax, Ay, Bx, By);
    if (x_inter <= x)
        return true;
    return false;
}

// Проверка, лежит ли точка левее прямой
bool is_left(double x, double y, double Ax, double Ay, double Bx, double By){
    double x_inter = get_x_inter(y, Ax, Ay, Bx, By);
    if (x_inter >= x)
        return true;
    return false;
}

// Вычисление коэффициента a_ij:
// Важно, что именно такой прямоугольник (симметрия относительно OY), и линия интеграла 
// может пересекать сторону треугольника только перпендикулярно оси OX.
double compute_a(int i, int j) {
    double x_left = A1 + (i - 0.5) * h1;
    double y_bottom = A2 + (j - 0.5) * h2;
    double y_top = A2 + (j + 0.5) * h2;

    double x_inter = abs(x_left);
    double y_inter = get_y_inter(x_inter, Ax, Ay, Bx, By); 

    if (y_top <= y_inter && y_bottom >= Cy)
        return 1.0;

    if (y_bottom >= y_inter || y_top <= Cy)
        return 1.0 / eps;

    double l = h2;
    if (y_bottom < Cy) {
        l -= Cy - y_bottom;
    }

    if (y_top > y_inter) {
        l -= y_top - y_inter;
    }
    
    double result = l/h2 + (1.0 - l/h2)/eps;
    
    return result;
}

// Вычисление коэффициента b_ij:
// Аналогично вычислению a_ij, только линия пересекает перпендикулярно OY, 
// а симметрией по OX не будем пользоваться, тк может быть двойной выход из области D.
double compute_b(int i, int j) {
    double x_left = A1 + (i - 0.5) * h1;
    double x_right = A1 + (i + 0.5) * h1;
    double y_bottom = A2 + (j - 0.5) * h2;

    double y_inter = y_bottom;
    double x_inter_left =  get_x_inter(y_inter, Cx, Cy, Bx, By);
    double x_inter_right =  get_x_inter(y_inter, Ax, Ay, Bx, By);

    if (x_right <= x_inter_right && x_left >= x_inter_left && y_bottom >= Cy)
        return 1.0;

    if (x_left >= x_inter_right || x_right <= x_inter_left || y_bottom < Cy)
        return 1.0 / eps;

    double l;
    if (x_left >= x_inter_left) {
        l = x_inter_right - x_left;
    } else if (x_right <= x_inter_right) {
        l = x_right - x_inter_left;
    } else {
        l = x_inter_right - x_inter_left; 
    }

    double result = l/h1 + (1.0 - l/h1)/eps;

    return result;
}

// Вычисление правой части F_ij
double compute_F(int i, int j) {
    // Граници прямоугольника, площадь вхождения
    // в D которого и надо определить 
    double x_left = A1 + (i - 0.5) * h1;
    double x_right = A1 + (i + 0.5) * h1;
    double y_bottom = A2 + (j - 0.5) * h2;
    double y_top = A2 + (j + 0.5) * h2;

    double s = h1*h2; 
    double out, x_inter, y_inter, x_inter_up, x_inter_down, x_inter_left, x_inter_right;

    // Прямоугольник вне области D
    if (is_left(x_right, y_bottom, Cx, Cy, Bx, By) 
        || is_right(x_left, y_bottom, Ax, Ay, Bx, By) || y_top <= Cy || y_bottom >= By
        || (y_bottom < Cy && is_left(x_right, y_top, Cx, Cy, Bx, By))
        || (y_bottom < Cy && is_right(x_left, y_top, Ax, Ay, Bx, By))) {
        
        return 0.0;
    }

    // Прямоугольник полностью в D
    if (is_left(x_right, y_top, Ax, Ay, Bx, By) 
        && is_right(x_left, y_top, Cx, Cy, Bx, By) && y_bottom >= Cy) {

        return 1.0;
    }

     // Нижний правый угол в области (как треугольник)
    if (is_right(x_right, y_bottom, Cx, Cy, Bx, By)
        &&  is_left(x_right, y_top, Cx, Cy, Bx, By) && y_bottom >= Cy && x_right <= Bx) {

        x_inter = get_x_inter(y_bottom, Cx, Cy, Bx, By);
        y_inter = get_y_inter(x_right, Cx, Cy, Bx, By);
        s = (x_right - x_inter) * (y_inter - y_bottom) / 2;
        return s /h1/h2;
    }

    // Нижний левый угол в области (как треугольник)
    if (is_left(x_left, y_bottom, Ax, Ay, Bx, By)
        &&  is_right(x_left, y_top, Ax, Ay, Bx, By) && y_bottom >= Cy && x_left >= Bx) {

        x_inter = get_x_inter(y_bottom, Ax, Ay, Bx, By);
        y_inter = get_y_inter(x_left, Ax, Ay, Bx, By);
        s = (x_inter - x_left) * (y_inter - y_bottom) / 2;
        return s /h1/h2;
    }

    // Подсчет площади в районе верхушки треугольника
    if (y_top > By && y_bottom < By) {
        // В области D находится треугольник
        if (is_right(x_right, y_bottom, Ax, Ay, Bx, By)
            &&  is_left(x_left, y_bottom, Cx, Cy, Bx, By)) {

            x_inter_left = get_x_inter(y_bottom, Cx, Cy, Bx, By);
            x_inter_right = get_x_inter(y_bottom, Ax, Ay, Bx, By);
            s = (By - y_bottom) * (x_inter_right - x_inter_left) / 2;
            return s / h1 / h2;
        }

        // В области D находится треугольник над прямоугольгиком
        if (is_left(x_right, y_bottom, Ax, Ay, Bx, By)
            &&  is_right(x_left, y_bottom, Cx, Cy, Bx, By)) {

            y_inter = get_y_inter(x_left, Cx, Cy, Bx, By);
            s = h1 * (y_inter - y_bottom);
            s += (By - y_inter) * h1 / 2;
            return s / h1 / h2;
        }

        // В области D находится треугольник над трапецией с острым нижним правым углом
        if (is_right(x_left, y_bottom, Cx, Cy, Bx, By)) {
            y_inter = get_y_inter(x_left, Cx, Cy, Bx, By);
            x_inter_up = get_x_inter(y_inter, Ax, Ay, Bx, By);
            x_inter_down = get_x_inter(y_bottom, Ax, Ay, Bx, By);
            s = (y_inter - y_bottom) * ((x_inter_up - x_left) + (x_inter_down - x_left)) / 2;
            s += (x_inter_up - x_left) * (By - y_inter) / 2;
            return s / h1 / h2;
        }

        // В области D находится треугольник над трапецией с острым нижним левым углом
        y_inter = get_y_inter(x_right, Ax, Ay, Bx, By);
        x_inter_up = get_x_inter(y_inter, Cx, Cy, Bx, By);
        x_inter_down = get_x_inter(y_bottom, Cx, Cy, Bx, By);
        s = (y_inter - y_bottom) * ((x_right - x_inter_up) + (x_right - x_inter_down)) / 2;
        s += (x_right - x_inter_up) * (By - y_inter) / 2;
        return s / h1 / h2;
    }

    // Вне области верхний левый угол (но как треугольник)
    if (is_left(x_left, y_top, Cx, Cy, Bx, By)
        && is_right(x_left, y_bottom, Cx, Cy, Bx, By) && y_top > Cy && y_top <= By) {
        
        x_inter = get_x_inter(y_top, Cx, Cy, Bx, By);
        y_inter = get_y_inter(x_left, Cx, Cy, Bx, By);
        out = (x_inter - x_left) * (y_top - y_inter) / 2;
        s -= out;
    }

    // Вне области верхний правый угол (но как треугольник)
    if (is_right(x_right, y_top, Ax, Ay, Bx, By)
        && is_left(x_right, y_bottom, Ax, Ay, Bx, By) && y_top > Cy && y_top <= By) {
        
        x_inter = get_x_inter(y_top, Ax, Ay, Bx, By);
        y_inter = get_y_inter(x_right, Ax, Ay, Bx, By);
        out = (x_right - x_inter) * (y_top - y_inter) / 2;
        s -= out;
    }

    // Вне области верхний левый угол (но как трапеция)
    if (is_left(x_left, y_bottom, Cx, Cy, Bx, By)
        && is_right(x_right, y_top, Cx, Cy, Bx, By) && y_top > Cy && y_top <= By) {
        
        x_inter_up = get_x_inter(y_top, Cx, Cy, Bx, By);
        x_inter_down = get_x_inter(y_bottom, Cx, Cy, Bx, By);
        out = ((x_inter_down - x_left) + (x_inter_up - x_left)) * h2 / 2;
        s -= out;
    }

    // Вне области верхний правый угол (но как трапеция)
    if (is_right(x_right, y_bottom, Ax, Ay, Bx, By) && y_top <= By
        && is_left(x_left, y_top, Ax, Ay, Bx, By) && y_top > Cy) {
        
        x_inter_up = get_x_inter(y_top, Ax, Ay, Bx, By);
        x_inter_down = get_x_inter(y_bottom, Ax, Ay, Bx, By);
        out = ((x_right - x_inter_down) + (x_right - x_inter_up)) * h2 / 2;
        s -= out;
    }

    // Вне области брюхо прямоугольника (но посередине треугольника)
    if (y_bottom < Cy && x_left >= Cx && x_right <= Ax) {

        out = (Cy - y_bottom) * (x_right - x_left);
        s -= out;
        return s /h1/h2;
    }

    // Вне области брюхо прямоугольника (но левым краем зашел в треугольник как трапеция)
    if (y_bottom < Cy && is_left(x_left, y_top, Ax, Ay, Bx, By)
        && is_right(x_right, y_top, Ax, Ay, Bx, By)) {

        x_inter_up = get_x_inter(y_top, Ax, Ay, Bx, By);
        x_inter_down = get_x_inter(Cy, Ax, Ay, Bx, By);

        s = (y_top - Cy) * ((x_inter_down - x_left) + (x_inter_up - x_left)) / 2;
        return s /h1/h2;
    }

    // Вне области брюхо прямоугольника (но правым краем зашел в треугольник как трапеция)
    if (y_bottom < Cy && is_right(x_right, y_top, Cx, Cy, Bx, By)
        && is_left(x_left, y_top, Cx, Cy, Bx, By)) {

        x_inter_up = get_x_inter(y_top, Cx, Cy, Bx, By);
        x_inter_down = get_x_inter(Cy, Cx, Cy, Bx, By);

        s = (y_top - Cy) * ((x_right - x_inter_down) + (x_right - x_inter_up)) / 2;
        return s /h1/h2;
    }

    return s/h1/h2; 
}

// Применение оператора A к сеточной функции
void applyOperator(vector<vector<double>>& Aw, const vector<vector<double>>& w, 
                   const vector<vector<double>>& a, const vector<vector<double>>& b) {
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            double term1 = a[i+1][j] * (w[i+1][j] - w[i][j]) - a[i][j] * (w[i][j] - w[i-1][j]);
            double term2 = b[i][j+1] * (w[i][j+1] - w[i][j]) - b[i][j] * (w[i][j] - w[i][j-1]);
            Aw[i][j] = -(term1 / (h1*h1) + term2 / (h2*h2));
        }
    }
}

// Решение системы Dz = r (диагональное предобуславливание)
void solveDiagonal( vector<vector<double>>& z, const vector<vector<double>>& r,
                   const vector<vector<double>>& a, const vector<vector<double>>& b) {
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            double diag = (a[i+1][j] + a[i][j]) / (h1 * h1) + (b[i][j+1] + b[i][j]) / (h2 * h2);
            z[i][j] = r[i][j] / diag;
        }
    }
}

// Скалярное произведение
double dot(const vector<vector<double>>& u, const vector<vector<double>>& v) {
    double result = 0.0;
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            result += u[i][j] * v[i][j];
        }
    }

    return result * h1 * h2;
}

// Норма разности сеточных функций
double normDifference(vector<vector<double>>& u_old, const vector<vector<double>>& u) {
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            u_old[i][j] -= u[i][j];
        }
    }
    return sqrt(dot(u_old, u_old));
}

void conjugateGradient(vector<vector<double>>& w, const vector<vector<double>>& a, 
                        const vector<vector<double>>& b, const vector<vector<double>>& F) {
    // Инициализация
    vector<vector<double>> r(M+1, vector<double>(N+1, 0.0));
    vector<vector<double>> z(M+1, vector<double>(N+1, 0.0));
    vector<vector<double>> p(M+1, vector<double>(N+1, 0.0));
    vector<vector<double>> Ap(M+1, vector<double>(N+1, 0.0));
    
    // Начальное приближение - нулевое
    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= N; j++) {
            w[i][j] = 0.0;
        }
    }
    
    // Вычисление начальной невязки
    applyOperator(Ap, w, a, b);
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            r[i][j] = F[i][j] - Ap[i][j];
        }
    }

    // Первая итерация (метод скорейшего спуска)
    solveDiagonal(z, r, a, b);

    // p1 = z0
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            p[i][j] = z[i][j];
        }
    }

    applyOperator(Ap, p, a, b);
    double alpha = dot(z, r) / dot(Ap, p);

    // Обновление решения
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            w[i][j] += alpha * p[i][j];
        }
    }

    // Основной итерационный процесс
    int max_iter = (M-1) * (N-1);
    vector<vector<double>> w_old(M+1, vector<double>(N+1, 0.0));
    
    for (int k = 1; k < max_iter; k++) {
        // Сохраняем предыдущее приближение
        w_old = w;
        double beta = 1.0 / dot(z, r);
        
        // Обновление невязки
        for (int i = 1; i < M; i++) {
            for (int j = 1; j < N; j++) {
                r[i][j] -= alpha * Ap[i][j];
            }
        }
        
        // Решение системы Dz = r
        solveDiagonal(z, r, a, b);
        
        // Вычисление коэффициента beta
        beta *= dot(z, r);
        
        // Новое направление спуска
        for (int i = 1; i < M; i++) {
            for (int j = 1; j < N; j++) {
                p[i][j] = z[i][j] + beta * p[i][j];
            }
        }
        
        // Вычисление Ap
        applyOperator(Ap, p, a, b);
        
        // Новый шаг спуска
        alpha = dot(z, r) / dot(Ap, p);
        
        // Обновление решения
        for (int i = 1; i < M; i++) {
            for (int j = 1; j < N; j++) {
                w[i][j] += alpha * p[i][j];
            }
        }
        
        // Проверка изменения решения
        double norm = normDifference(w_old, w);
        if (norm < delta) {
            cout << "Сходимость достигнута на итерации " << k << endl;
            break;
        }
        
        // if (k % 10 == 0) {
        //     cout << "Итерация " << k << ", невязка: " << sqrt(dot(r, r)) << endl;
        // }

        if (k % 100 == 0) {
            cout << "Итерация " << k << ", норма разности: " << norm << endl;
        }
    }

}

int main() {
    cout << "Решение задачи Дирихле для уравнения Пуассона в треугольнике" << endl;
    cout << "Сетка: " << M << " x " << N << endl;
    cout << "Шаги: h1 = " << h1 << ", h2 = " << h2 << endl;
    cout << "Параметр epsilon = " << eps << endl;
    cout << "Параметр delta = " << delta << endl << endl;

    // Инициализация массивов
    vector<vector<double>> a(M+1, vector<double>(N+1, 0.0));
    vector<vector<double>> b(M+1, vector<double>(N+1, 0.0));
    vector<vector<double>> F(M, vector<double>(N, 0.0));
    vector<vector<double>> w(M+1, vector<double>(N+1, 0.0));

    clock_t start_time = clock();
    
    cout << "Начало инициализации..." << endl;

    // Инициализация коэффициентов
    for (int i = 1; i <= M; i++) {
        for (int j = 1; j <= N; j++) {
            a[i][j] = compute_a(i, j);
            b[i][j] = compute_b(i, j);
        }
    }

    // Инициализация правой части
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            F[i][j] = compute_F(i, j);
        }
    }

    double init_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    cout << "Время инициализации: " << init_time << " секунд" << endl;

    
    clock_t solve_start = clock();

    // Решение системы методом сопряженных градиентов
    conjugateGradient(w, a, b, F);

    double solve_time = (double)(clock() - solve_start) / CLOCKS_PER_SEC;
    cout << "Решение СЛАУ завершено за " << solve_time << " секунд" << endl << endl;

    // Общее время выполнения
    double total_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;

    cout << "Время решения СЛАУ: " << solve_time << " секунд" << endl;
    cout << "Общее время выполнения: " << total_time << " секунд" << endl;
    
    ofstream outfile("solution.txt");
    outfile << "x y w" << endl;
    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= N; j++) {
            double x = A1 + i * h1;
            double y = A2 + j * h2;
            outfile << x << " " << y << " " << w[i][j] << endl;
        }
    }
    outfile.close();
    
    cout << "Результаты записаны в solution.txt" << endl;
    

    return 0;
}