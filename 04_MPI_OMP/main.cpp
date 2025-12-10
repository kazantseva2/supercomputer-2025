#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime>
#include <omp.h>
#include <mpi.h>

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
const int M = 800, N = 1200;
const double h1 = (B1 - A1) / M;
const double h2 = (B2 - A2) / N;

const double h = max(h1, h2);
const double eps = h*h;    // малый параметр метода фиктивных областей
const double delta = 1e-8;      // точность итерационного метода

struct  Domain {
    int start_i, end_i;
    int start_j, end_j;
    int size_M, size_N;

    int left, right;
    int bottom, top;
    
    int rank;
    MPI_Comm comm_cart;

    Domain() {
        int size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // Создаем декартову топологию 2×2
        int dims[2] = {0, 0};
        MPI_Dims_create(size, 2, dims);  // Автоматическое разбиение

        // Меняем местами, так как узлов больше по Y
        if (N > M) {
            int t = dims[0];
            dims[0] = dims[1];
            dims[1] = t;
        }

        int periods[2] = {0, 0};  // Периодические границы
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_cart);

        // Узнаем свои координаты
        int coords[2];
        MPI_Cart_coords(comm_cart, rank, 2, coords);
        
        // Находим соседей
        MPI_Cart_shift(comm_cart, 0, 1, &left, &right);   // По X
        MPI_Cart_shift(comm_cart, 1, 1, &bottom, &top);   // По Y

        if (rank == 0) {
            cout << "Разбиение:" << dims[0] << "*" << dims[1] << endl;
        }

        int m = (M - 1) / dims[0];
        int n = (N - 1) / dims[1];
        int m_ost = (M - 1) % dims[0];
        int n_ost = (N - 1) % dims[1];

        start_i = 1 + coords[0] * m + min(coords[0], m_ost);
        start_j = 1 + coords[1] * n + min(coords[1], n_ost);

        end_i = start_i + m - 1 + (coords[0] < m_ost ? 1 : 0);
        end_j = start_j + n - 1 + (coords[1] < n_ost ? 1 : 0);

        size_M = end_i - start_i + 1;
        size_N = end_j - start_j + 1;

        // cout << "Процесс " << rank << "(" << coords[0] << "," << coords[1] 
        // << "), start_i=" << start_i << ", end_i=" << end_i
        // << ", start_j=" << start_j << ", end_j=" << end_j 
        // << ", size_M=" << size_M << ", size_N=" << size_N << endl;
    }
};

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

    // Вне области брюво прямоугольника (но левым краем зашел в треугольник как трапеция)
    if (y_bottom < Cy && is_left(x_left, y_top, Ax, Ay, Bx, By)
        && is_right(x_right, y_top, Ax, Ay, Bx, By)) {

        x_inter_up = get_x_inter(y_top, Ax, Ay, Bx, By);
        x_inter_down = get_x_inter(Cy, Ax, Ay, Bx, By);

        s = (y_top - Cy) * ((x_inter_down - x_left) + (x_inter_up - x_left)) / 2;
        return s /h1/h2;
    }

    // Вне области брюво прямоугольника (но правым краем зашел в треугольник как трапеция)
    if (y_bottom < Cy && is_right(x_right, y_top, Cx, Cy, Bx, By)
        && is_left(x_left, y_top, Cx, Cy, Bx, By)) {

        x_inter_up = get_x_inter(y_top, Cx, Cy, Bx, By);
        x_inter_down = get_x_inter(Cy, Cx, Cy, Bx, By);

        s = (y_top - Cy) * ((x_right - x_inter_down) + (x_right - x_inter_up)) / 2;
        return s /h1/h2;
    }

    return s/h1/h2; 
}

// Функция для инициализации коэффициентов a и b с вычислением граничных значений
void initializeCoefficients(vector<vector<double>>& a, vector<vector<double>>& b, Domain d) {
    #pragma omp parallel for collapse(2)
    for (int i_loc = 0; i_loc < d.size_M + 2; i_loc++) {
        for (int j_loc = 0; j_loc < d.size_N + 2; j_loc++) {
            // Вычисляем глобальные индексы
            int i_global = d.start_i - 1 + i_loc;  // -1 потому что i_loc=0 это левый дополнительный край
            int j_global = d.start_j - 1 + j_loc;  // аналогично

            // Нижняя строка и левый столбец не заполняются по алгоритму из задания
            if (i_global != 0 && j_global != 0) {
                a[i_loc][j_loc] = compute_a(i_global, j_global);
                b[i_loc][j_loc] = compute_b(i_global, j_global);
            }
        }
    }
}

// Функция для инициализации правой части F (только внутренняя область)
void initializeF(vector<vector<double>>& F, Domain d) {
    #pragma omp parallel for collapse(2)
    for (int i_loc = 1; i_loc <= d.size_M ; i_loc++) {
        for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
            int i_global = d.start_i - 1 + i_loc;
            int j_global = d.start_j - 1 + j_loc;

            F[i_loc][j_loc] = compute_F(i_global, j_global);
        }
    }
}

void applyOperator(vector<vector<double>>& Aw, vector<vector<double>>& w, 
                   const vector<vector<double>>& a, const vector<vector<double>>& b,
                   Domain d) {

    // Создаем буферы для граничных значений
    vector<double> send_left(d.size_N, 0.0), send_right(d.size_N, 0.0);
    vector<double> send_bottom(d.size_M, 0.0), send_top(d.size_M, 0.0);
    vector<double> recv_left(d.size_N, 0.0), recv_right(d.size_N, 0.0);
    vector<double> recv_bottom(d.size_M, 0.0), recv_top(d.size_M, 0.0);
    
    // Заполняем буферы для отправки из ВНУТРЕННИХ границ
    for (int j = 0; j < d.size_N; j++) {
        send_left[j] = w[1][j+1];        // Внутренняя левая граница 
        send_right[j] = w[d.size_M][j+1]; // Внутренняя правая граница
    }
    
    for (int i = 0; i < d.size_M; i++) {
        send_bottom[i] = w[i+1][1];      // Внутренняя нижняя граница
        send_top[i] = w[i+1][d.size_N];  // Внутренняя верхняя граница
    }

    // Неблокирующие операции для обмена границами
    MPI_Request requests[8] = {MPI_REQUEST_NULL};
    int request_count = 0;

    // ЛЕВЫЙ СОСЕД
    if (d.left != MPI_PROC_NULL) {
        MPI_Isend(send_left.data(), d.size_N, MPI_DOUBLE, d.left, 0, d.comm_cart, &requests[request_count++]);
        MPI_Irecv(recv_left.data(), d.size_N, MPI_DOUBLE, d.left, 0, d.comm_cart, &requests[request_count++]);
    } 

    // ПРАВЫЙ СОСЕД
    if (d.right != MPI_PROC_NULL) {
        MPI_Isend(send_right.data(), d.size_N, MPI_DOUBLE, d.right, 0, d.comm_cart, &requests[request_count++]);
        MPI_Irecv(recv_right.data(), d.size_N, MPI_DOUBLE, d.right, 0, d.comm_cart, &requests[request_count++]);
    }

    // НИЖНИЙ СОСЕД
    if (d.bottom != MPI_PROC_NULL) {
        MPI_Isend(send_bottom.data(), d.size_M, MPI_DOUBLE, d.bottom, 0, d.comm_cart, &requests[request_count++]);
        MPI_Irecv(recv_bottom.data(), d.size_M, MPI_DOUBLE, d.bottom, 0, d.comm_cart, &requests[request_count++]);
    } 

    // ВЕРХНИЙ СОСЕД
    if (d.top != MPI_PROC_NULL) {
        MPI_Isend(send_top.data(), d.size_M, MPI_DOUBLE, d.top, 0, d.comm_cart, &requests[request_count++]);
        MPI_Irecv(recv_top.data(), d.size_M, MPI_DOUBLE, d.top, 0, d.comm_cart, &requests[request_count++]);
    }

    // Ждем завершения всех коммуникаций 
    MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);

    // Заполняем полученными значениями
    for (int j = 0; j < d.size_N; j++) {
        w[0][j+1] = recv_left[j];        // Левый
        w[d.size_M+1][j+1] = recv_right[j]; // Правый
    }
    
    for (int i = 0; i < d.size_M; i++) {
        w[i+1][0] = recv_bottom[i];      // Нижний
        w[i+1][d.size_N+1] = recv_top[i]; // Верхний
    }

    // Теперь применяем оператор ко ВСЕМ внутренним точкам 
    #pragma omp parallel for collapse(2)
    for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
        for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
            
            // Используем разностную схему
            double term1 = a[i_loc+1][j_loc] * (w[i_loc+1][j_loc] - w[i_loc][j_loc]) 
                         - a[i_loc][j_loc] * (w[i_loc][j_loc] - w[i_loc-1][j_loc]);
            
            double term2 = b[i_loc][j_loc+1] * (w[i_loc][j_loc+1] - w[i_loc][j_loc]) 
                         - b[i_loc][j_loc] * (w[i_loc][j_loc] - w[i_loc][j_loc-1]);
            
            Aw[i_loc][j_loc] = -(term1 / (h1*h1) + term2 / (h2*h2));
        }
    }
}

// Решение системы Dz = r (диагональное предобуславливание)
void solveDiagonal(vector<vector<double>>& z, const vector<vector<double>>& r,
                   const vector<vector<double>>& a, const vector<vector<double>>& b,
                   Domain d) {
    #pragma omp parallel for collapse(2)
    for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
        for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
            double diag = (a[i_loc+1][j_loc] + a[i_loc][j_loc]) / (h1 * h1) 
                        + (b[i_loc][j_loc+1] + b[i_loc][j_loc]) / (h2 * h2);
            z[i_loc][j_loc] = r[i_loc][j_loc] / diag;
        }
    }
}

// Скалярное произведение с MPI редукцией
double dot(const vector<vector<double>>& u, const vector<vector<double>>& v, Domain d) {
    double local_result = 0.0;
    
    // Суммируем только по внутренней области процесса
    #pragma omp parallel for collapse(2) reduction(+:local_result)
    for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
        for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
            local_result += u[i_loc][j_loc] * v[i_loc][j_loc];
        }
    }
    local_result *= h1 * h2; 
    
    // Суммируем по всем процессам
    double global_result = 0.0;
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return global_result;
}

// Норма разности сеточных функций
double normDifference(vector<vector<double>>& u_old, const vector<vector<double>>& u, Domain d) {
    double local_norm = 0.0;
    
    #pragma omp parallel for collapse(2) reduction(+:local_norm)
    for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
        for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
            double diff = u_old[i_loc][j_loc] - u[i_loc][j_loc];
            local_norm += diff * diff;
        }
    }
    
    local_norm *= h1 * h2;
    
    // Суммируем по всем процессам
    double global_norm = 0.0;
    MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return sqrt(global_norm);
}

void conjugateGradient(vector<vector<double>>& w, const vector<vector<double>>& a, 
                    const vector<vector<double>>& b, const vector<vector<double>>& F,
                    Domain d) {

    // Инициализация рабочих массивов
    vector<vector<double>> r(d.size_M + 1, vector<double>(d.size_N + 1, 0.0));
    vector<vector<double>> z(d.size_M + 1, vector<double>(d.size_N + 1, 0.0));
    vector<vector<double>> p(d.size_M + 2, vector<double>(d.size_N + 2, 0.0));
    vector<vector<double>> Ap(d.size_M + 1, vector<double>(d.size_N + 1, 0.0));
    vector<vector<double>> w_old(d.size_M + 2, vector<double>(d.size_N + 2, 0.0));

    // Вычисление начальной невязки
    applyOperator(Ap, w, a, b, d);
    
    #pragma omp parallel for collapse(2)
    for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
        for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
            r[i_loc][j_loc] = F[i_loc][j_loc] - Ap[i_loc][j_loc];
        }
    }

    // Первая итерация (метод скорейшего спуска)
    solveDiagonal(z, r, a, b, d);

    // p1 = z0
    #pragma omp parallel for collapse(2)
    for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
        for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
            p[i_loc][j_loc] = z[i_loc][j_loc];
        }
    }

    applyOperator(Ap, p, a, b, d);
    double alpha = dot(z, r, d) / dot(Ap, p, d);

    // Обновление решения
    #pragma omp parallel for collapse(2)
    for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
        for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
            w[i_loc][j_loc] += alpha * p[i_loc][j_loc];
        }
    }

    // Основной итерационный процесс
    int max_iter = d.size_M * d.size_N;  // Локальный размер

    for (int k = 1; k < max_iter; k++) {
        // Сохраняем предыдущее приближение
        w_old = w;

        double beta = 1.0 / dot(z, r, d);

        // Обновление невязки
        #pragma omp parallel for collapse(2)
        for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
            for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
                r[i_loc][j_loc] -= alpha * Ap[i_loc][j_loc];
            }
        }

        // Решение системы Dz = r
        solveDiagonal(z, r, a, b, d);

        // Вычисление коэффициента beta
        beta *= dot(z, r, d);

        // Новое направление спуска
        #pragma omp parallel for collapse(2)
        for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
            for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
                p[i_loc][j_loc] = z[i_loc][j_loc] + beta * p[i_loc][j_loc];
            }
        }

        // Вычисление Ap
        applyOperator(Ap, p, a, b, d);
        
        // Новый шаг спуска
        alpha = dot(z, r, d) / dot(Ap, p, d);

        // Обновление решения
        #pragma omp parallel for collapse(2)
        for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
            for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
                w[i_loc][j_loc] += alpha * p[i_loc][j_loc];
            }
        }

        // Проверка изменения решения (синхронизировано по всем процессам)
        double current_norm = normDifference(w_old, w, d);
        
        if (d.rank == 0 && k % 100 == 0) {
            cout << "Итерация " << k << ", норма изменения: " << current_norm << endl;
        }
        
        if (current_norm < delta) {
            if (d.rank == 0) {
                cout << "Сходимость достигнута на итерации " << k << endl;
            }
            break;
        }
    }
}

void save(vector<vector<double>>& w, Domain d) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (d.rank == 0) {
        ofstream outfile("solution.txt");
        outfile << "x y w" << endl;
        outfile.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    // Процессы пишут по очереди
    for (int proc = 0; proc < size; proc++) {
        if (d.rank == proc) {
            ofstream outfile("solution.txt", ios::app);
            for (int i_loc = 1; i_loc <= d.size_M; i_loc++) {
                for (int j_loc = 1; j_loc <= d.size_N; j_loc++) {
                    int i_global = d.start_i - 1 + i_loc;
                    int j_global = d.start_j - 1 + j_loc;
                    double x = A1 + i_global * h1;
                    double y = A2 + j_global * h2;
                    outfile << x << " " << y << " " << w[i_loc][j_loc] << endl;
                }
            }
            outfile.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    if (d.rank == 0) {
        cout << "Результаты записаны в solution.txt" << endl;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Установка числа потоков OpenMP
    int num_threads = 1;
    if (argc > 1) {
        num_threads = atoi(argv[1]);
    }
    omp_set_num_threads(num_threads);

    double start_time, init_time, solve_time, total_time;
    
    if (rank == 0) {
        cout << "Решение задачи Дирихле для уравнения Пуассона в треугольнике" << endl;
        cout << "Сетка: " << M << " x " << N << endl;
        cout << "Шаги: h1 = " << h1 << ", h2 = " << h2 << endl;
        cout << "Параметр epsilon = " << eps << endl;
        cout << "Параметр delta = " << delta << endl << endl;
        cout << "Количество MPI процессов: " << size << endl;
        cout << "Количество OpenMP потоков на процесс: " << num_threads << endl;
        cout << "Общее число потоков: " << size * num_threads << endl << endl;
    }
    
    start_time = MPI_Wtime();

    Domain d = Domain();

    // Инициализация массивов
    vector<vector<double>> a(d.size_M + 2, vector<double>(d.size_N + 2, 0.0));
    vector<vector<double>> b(d.size_M + 2, vector<double>(d.size_N + 2, 0.0));
    vector<vector<double>> F(d.size_M + 2, vector<double>(d.size_N + 2, 0.0));
    vector<vector<double>> w(d.size_M + 2, vector<double>(d.size_N + 2, 0.0));
    
    if (rank == 0) {
        cout << "Начало инициализации..." << endl;
    }

    initializeCoefficients(a, b, d);
    initializeF(F, d);

    init_time = MPI_Wtime() - start_time;
    
    if (rank == 0) {
        cout << "Решаем СЛАУ..." << endl;
    }

    // Синхронизация
    MPI_Barrier(MPI_COMM_WORLD);
    solve_time = MPI_Wtime();

    // Решение системы методом сопряженных градиентов
    conjugateGradient(w, a, b, F, d);

    MPI_Barrier(MPI_COMM_WORLD);

    solve_time = MPI_Wtime() - solve_time;
    total_time = MPI_Wtime() - start_time;
    
    // Сбор результатов
    if (rank == 0) {
        cout << endl << "Время инициализации: " << init_time << " секунд" << endl;
        cout << "Время решения СЛАУ: " << solve_time << " секунд" << endl;
        cout << "Общее время выполнения: " << total_time << " секунд" << endl << endl;
    }

    //save(w, d);
        
    MPI_Finalize();
    return 0;
}