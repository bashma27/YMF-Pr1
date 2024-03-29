﻿#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <set>
#include <iomanip>
#include <math.h>
#define EPS 1e-16 //Для сравнения double

using namespace std;

#pragma region Объявление глобальных переменных
int NUM__NODES_IN_EDGE_X, NUM__NODES_IN_EDGE_Y, NUM_NODES; //число узлов в ребре по x/y, число узлов
int Nx, Ny, L; // Nx - число вертикальных границ, Ny - число горизонтальных границ, L - число прямоугольников
double LAMBDA = 1; // коеф. лямбда
double GAMMA = 2; // коеф. гамма
int NUM_BOUND_COND_1, NUM_BOUND_COND_2; //Количество рёбер с первыми и вторыми краевыми условиями
int m; // расстояние между диагоналями
vector<vector<double>> al, au;
vector<double> Xw, Yw; // границы области
vector<vector<int>> W; // прямоугольники
vector<double> di, q, b;
vector<pair<double, double>> nodes; // узлы сетки, заданные координатами 
vector<vector<int>> edge; // первый индекс номер краевого, потом индексы
vector<vector<pair<int, double>>> koef; // Коэффициеты n, q  для каждой области в виде пары

#pragma endregion

#pragma region Работа с сеткой 
int InputBorders() {
    ifstream File_Grid("grid.txt");
    if (!File_Grid.is_open()) return 1;
    File_Grid >> Nx;
    Xw.resize(Nx);
    for (int i = 0; i < Nx; i++) {
        File_Grid >> Xw[i];
    }
    File_Grid >> Ny;
    Yw.resize(Ny);
    for (int i = 0; i < Ny; i++) {
        File_Grid >> Yw[i];
    }
    File_Grid >> L;
    W.resize(L);
    for (int i = 0; i < L; i++)
        W[i].resize(4);
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < 4; j++) {
            File_Grid >> W[i][j];
        }
    }
    File_Grid.close();
    return 0;
}

int CalcSum_ni(int index, int start, int stop) { // Посчитать сумму ni по порядку от start до end
    int value = 0;
    for (int k = start; k < stop; k++) {
        value += koef[index][k].first;
    }
    return value;
}

//ФУНКЦИЯ ОПРЕДЕЛЯЮЩАЯ ФИКТИВНЫЙ УЗЕЛ
bool IsFictitious(int num) {
    double x = nodes[num].first, y = nodes[num].second;
    int x1, x2, y1, y2;
    for (int k = 0; k < L; k++) {
        x1 = W[k][0]; y1 = W[k][2]; //берём индексы координат прямоугольника
        x2 = W[k][1]; y2 = W[k][3];
        if (Xw[x1 - 1] <= x && Xw[x2 - 1] >= x && Yw[y1 - 1] <= y && Yw[y2 - 1] >= y) // проверка на то что узел находится в прямоугольнике
            return false;
    }
    return true; // если по результату прохода узел не находится не в одном прямоугольнике то фиктивный узел.
}

int Input_koef() {
    ifstream File_Koef("koef.txt");
    if (!File_Koef.is_open()) return 1;
    koef.resize(2);
    koef[0].resize(Nx - 1);
    koef[1].resize(Ny - 1);
    int ni;
    double qi;
    for (int i = 0; i < Nx - 1; i++) {
        File_Koef >> ni >> qi;
        koef[0][i] = { ni, qi };
    }
    for (int i = 0; i < Ny - 1; i++) {
        File_Koef >> ni >> qi;
        koef[1][i] = { ni, qi };
    }
    File_Koef.close();
    NUM__NODES_IN_EDGE_X = CalcSum_ni(0, 0, Nx - 1) + 1;
    NUM__NODES_IN_EDGE_Y = CalcSum_ni(1, 0, Ny - 1) + 1;
    NUM_NODES = NUM__NODES_IN_EDGE_X * NUM__NODES_IN_EDGE_Y;
    nodes.resize(NUM_NODES);
    return 0;
}

bool InVector(int value, vector<int>& temp) { //Для записи в более коротком формате.
    for (auto var : temp) {
        if (var == value)
            return false;
    }
    return true;
}
#pragma endregion

#pragma region Задание функций

double f(double x, double y) {
    return 2 * (x * x * x * x + y * y * y * y) - 12 * (x * x + y * y);
}

double theta(double x, double y) { // краевое условие второго рода
    return - 8;
}

double u_g(double x, double y) { // краевое условие первого рода
    return x * x * x * x + y * y * y * y;
}

#pragma endregion

#pragma region Вывод вектора решений
void Output() {
    ofstream f_result;
    f_result.imbue(locale("Russian"));
    f_result.open("f_result.txt");
    for (int i = 0; i < NUM_NODES; i++) {
        f_result << setprecision(15) << q[i] << endl;
    }
    f_result.close();
}
#pragma endregion

#pragma region Генерация сетки
int Create_Grid() {
    InputBorders(); // Загрузка границ области и прямоугольников
    Input_koef(); // Загрузка коэфициентов и объявление некторых констант
    int nx = 0, ny = 0, j = 0, i = 0; //nx и ny для подсчёта кол-ва разбиений до нужной точки
    int it = 0; // текущий номер последнего элемента в nodes
    double qx = 1.0, qy = 1.0, hy = 0.0, hx = 0.0, y = 0.0;
    for (i = 0; i < Ny - 1; i++) { // перебираем все коэфициенты по y
        ny = koef[1][i].first;
        qy = koef[1][i].second;
        if (abs(qy - 1.0) < EPS)
            hy = (Yw[i + 1] - Yw[i]) / ny; // вычисления шага для равномерной сетки
        else
            hy = (Yw[i + 1] - Yw[i]) * (qy - 1) / (pow(qy, ny) - 1); //вычисление h0 для неравномерной сетки
        for (int l = 0; l < ny; l++, it++) { // цикл по каждому разбиению по y
            if (abs(qy - 1.0) < EPS) // подсчёт текущего y
                y = Yw[i] + hy * l;
            else
                y = Yw[i] + hy * (pow(qy, l) - 1) / (qy - 1);
            for (j = 0; j < Nx - 1; j++) { // перебор всех коэфициентов по x
                nx = koef[0][j].first;
                qx = koef[0][j].second;
                if (abs(qx - 1.0) < EPS)
                    hx = (Xw[j + 1] - Xw[j]) / nx; // вычисление шага для равномерной сетки
                else
                    hx = (Xw[j + 1] - Xw[j]) * (qx - 1) / (pow(qx, nx) - 1); // вычисление h0 для неравномерной сетки
                for (int k = 0; k < nx; k++, it++) {
                    if (abs(qx - 1.0) < EPS) // считаем текущий x и кладём всё в nodes
                        nodes[it] = { Xw[j] + k * hx, y };
                    else
                        nodes[it] = { Xw[j] + hx * (pow(qx, k) - 1) / (qx - 1), y };
                }
            }
            nodes[it] = { Xw[j], y }; // Кладётся конец для x, иначе бы элемент бы повторялся
        }

    }
    //---------------------------------------------------------
    // Дальше идёт тоже самое для последней строчки
    //---------------------------------------------------------
    y = Yw[i];
    for (j = 0; j < Nx - 1; j++) {
        nx = koef[0][j].first;
        qx = koef[0][j].second;
        if (abs(qx - 1.0) < EPS)
            hx = (Xw[j + 1] - Xw[j]) / nx;
        else
            hx = (Xw[j + 1] - Xw[j]) * (qx - 1) / (pow(qx, nx) - 1);
        for (int k = 0; k < nx; k++, it++) {
            if (abs(qx - 1.0) < EPS)
                nodes[it] = { Xw[j] + k * hx, y };
            else
                nodes[it] = { Xw[j] + hx * (pow(qx, k) - 1) / (qx - 1), y };
        }
    }
    nodes[it] = { Xw[j], y };
    //---------------------------------------------------------------

    // Заполнение edge
    /*
        Пример ввода edge.txt (Сначало 1 краевые, потом вторые краевые)
        6 6
        1 4 1 1
        4 4 1 2
        1 1 1 6
        4 4 5 6
        ...
        1 4 6 6
    */
    ifstream File_Edge("edge.txt");
    if (!File_Edge.is_open()) return 1;
    int x1 = 0, x2 = 0, y1 = 0, y2 = 0, ij = 0, side = 0;
    File_Edge >> NUM_BOUND_COND_1 >> NUM_BOUND_COND_2; //считывание кол-во ребер каждого краевого условия
    edge.resize(2);
    vector<int> corner(NUM_BOUND_COND_1 + NUM_BOUND_COND_2, -1);
    it = 0;
    for (int l = 0; l < NUM_BOUND_COND_1; l++) { // Сначало обрабатываем 1 краевые
        File_Edge >> x1 >> x2 >> y1 >> y2;
        // Находим номер первого узла ребра (для наглядности разделил на две операции по каждому индексу)
        i = CalcSum_ni(0, 0, x1 - 1);
        j = CalcSum_ni(1, 0, y1 - 1) * NUM__NODES_IN_EDGE_X;
        ij = i + j; // номер первого ребра
        if (x2 - x1 > 0) { // Ребро горизонтальное
            nx = CalcSum_ni(0, x1 - 1, x2 - 1); // Кол-во разбиений на этом ребре по x
            if (InVector(ij, corner)) {
                edge[0].push_back(ij);
                corner[it] = ij;
                it++;
            }
            for (int k = 1; k < nx; k++) //Тут в связи с структурой данных забиваем на повторы (set сам уберет повторы)
                edge[0].push_back(ij + k);
            if (InVector(ij + nx, corner)) {
                edge[0].push_back(ij + nx);
                corner[it] = ij + nx;
                it++;
            }
        }
        else {
            ny = CalcSum_ni(1, y1 - 1, y2 - 1); //Кол-во разбиений на этом ребре по y          
            if (InVector(ij, corner)) {
                edge[0].push_back(ij);
                corner[it] = ij;
                it++;
            }
            for (int k = 1; k < ny; k++) // Так же забиваем на повторы
                edge[0].push_back(ij + k * NUM__NODES_IN_EDGE_X);
            if (InVector(ij + ny * NUM__NODES_IN_EDGE_X, corner)) {
                edge[0].push_back(ij + ny * NUM__NODES_IN_EDGE_X);
                corner[it] = ij + ny * NUM__NODES_IN_EDGE_X;
                it++;
            }
        }
    }
    for (int l = 0; l < NUM_BOUND_COND_2; l++) { //Обрабатываем вторые краевыеы
        File_Edge >> x1 >> x2 >> y1 >> y2;
        //Находим номер первого узла ребра (для наглядности разделил на две операции по каждому индексу)
        i = CalcSum_ni(0, 0, x1 - 1);
        j = CalcSum_ni(1, 0, y1 - 1) * NUM__NODES_IN_EDGE_X;
        ij = i + j; //номер первого узла ребра
        if (x2 - x1 > 0) { // Ребро горизонтальное
            nx = CalcSum_ni(0, x1 - 1, x2 - 1); // Кол-во разбиений на этом ребре по x           
            if (InVector(ij, corner)) {
                edge[1].push_back(ij);
                corner[it] = ij;
                it++;
            }
            for (int k = 1; k < nx; k++) // Эти не совпадают
                edge[1].push_back(ij + k);
            if (InVector(ij + nx, corner)) {
                edge[1].push_back(ij + nx);
                corner[it] = ij + nx;
                it++;
            }
        }
        else {
            ny = CalcSum_ni(1, y1 - 1, y2 - 1); // Кол-во разбиений на этом ребре по y
            if (InVector(ij, corner)) {
                edge[1].push_back(ij);
                corner[it] = ij;
                it++;
            }
            for (int k = 1; k < ny; k++) // Эти не совпадают
                edge[1].push_back(ij + k * NUM__NODES_IN_EDGE_X);
            if (InVector(ij + ny * NUM__NODES_IN_EDGE_X, corner)) {
                edge[1].push_back(ij + ny * NUM__NODES_IN_EDGE_X);
                corner[it] = ij + ny * NUM__NODES_IN_EDGE_X;
                it++;
            }
        }
    }
    File_Edge.close();
    corner.clear();
}
#pragma endregion

#pragma region Работа с краевыми условиями и фиктивными узлами

void ConsiderBoundConditFirstType() { // учет краевых условий первого типа
   
    for (int i = 0; i < edge[0].size(); i++) {
        int j = edge[0][i]; 
        b[j] = u_g(nodes[j].first, nodes[j].second);
        di[j] = double(1);
    }

}

void ConsiderBoundConditSecType() { // учет краевых условий второго типа
    for (int i = 0; i < edge[1].size(); i++) {
        int j = edge[1][i];
        int up_point = j + NUM__NODES_IN_EDGE_X;
        au[0][up_point] = -LAMBDA / (nodes[up_point].second - nodes[j].second);
        di[j] = LAMBDA / (nodes[up_point].second - nodes[j].second);
        b[j] = theta(nodes[j].first, nodes[j].second);
    }
}

void ConsiderFictitiousNodes() { // учет фиктивных узлов
    for (int i = 0; i < NUM_NODES; i++) {
        if (IsFictitious(i)) {
            b[i] = 0;
            di[i] = double(1);
        }
    }
}
#pragma endregion

#pragma region Построение глобальной матрицы и вектора

void BuildMatrA() {
    al.assign(2, vector<double> (NUM_NODES));
    au.assign(2, vector<double> (NUM_NODES));
    di.resize(NUM_NODES);
    for (int i = 0; i < NUM_NODES; i++) {
        if (!InVector(i, edge[0]) || !InVector(i, edge[1])) continue;
        if (IsFictitious(i)) continue;

        int left_point, right_point, down_point, up_point;
        double h_x_prev, h_x_curr, h_y_prev, h_y_curr; // предыдущий / текущий
        down_point = i - NUM__NODES_IN_EDGE_X;
        up_point = i + NUM__NODES_IN_EDGE_X;
        left_point = i - 1;
        right_point = i + 1;
        h_x_prev = nodes[i].first - nodes[left_point].first;
        h_x_curr = nodes[right_point].first - nodes[i].first;
        h_y_prev = nodes[i].second - nodes[down_point].second;
        h_y_curr = nodes[up_point].second - nodes[i].second;
        al[0][i] = -2 * LAMBDA / (h_y_prev * (h_y_curr + h_y_prev));
        al[1][i] = -2 * LAMBDA / (h_x_prev * (h_x_curr + h_x_prev));
        au[0][i + NUM__NODES_IN_EDGE_X] = -2 * LAMBDA / (h_y_curr * (h_y_curr + h_y_prev));
        au[1][i + 1] = -2 * LAMBDA / (h_x_curr * (h_x_curr + h_x_prev));
        di[i] = 2 * LAMBDA * (1 / (h_x_curr * h_x_prev) + 1 / (h_y_curr * h_y_prev)) + GAMMA;
        
    }
}

void BuildVecB() {
    b.resize(NUM_NODES);
    for (int i = 0; i < NUM_NODES; i++) {
        if (!InVector(i, edge[0]) || !InVector(i, edge[1])) continue;
        if (IsFictitious(i)) continue;
        b[i] = f(nodes[i].first, nodes[i].second);
    }
}

#pragma endregion

#pragma region Решение СЛАУ(Гаусс-Зейдель)

double norma() {
    double res = 0;
    for (int i = 0; i < NUM_NODES; i++) {
        res += b[i] * b[i];
    }
    return sqrt(res);
}

double iter(int i) {
    double temp = 0;
    temp += di[i] * q[i];

    if (i > 0) {
        temp += al[1][i] * q[i - 1];
    }
    if (i > m + 1) {
        temp += al[0][i] * q[i - m - 2];
    }

    if (i < NUM_NODES - 1) {
        temp += au[1][i + 1] * q[i + 1];
    }
    if (i < NUM_NODES - m - 2) {
        temp += au[0][i + m + 2] * q[i + m + 2];
    }
    temp = b[i] - temp;

    return temp;
}

void GaussZaid(double w, double eps, int max_iter) {
    int counter = 1;
    double resid = 1;
    double norm = norma();
    for (; counter < max_iter && resid > eps; counter++) {
        resid = 0;
        for (int i = 0; i < NUM_NODES; i++) {
            double temp = iter(i);
            q[i] = q[i] + w * temp / di[i];
            resid += temp * temp;
        }
        resid = sqrt(resid) / norm;
        cout << setprecision(15) << resid << " " << counter << endl;
    }
}
#pragma endregion

#pragma region Тестирование
void Test() {
    vector<double> q_u(NUM_NODES, 0);
    for (int i = 0; i < NUM_NODES; i++) {
        q_u[i] = u_g(nodes[i].first, nodes[i].second);
    }
    double norm_vec_err = 0, norm_vec_q_u = 0; // норма вектора погрешности и q_u
    for (int i = 0; i < NUM_NODES; i++) {
        if (IsFictitious(i)) continue;
        norm_vec_err += (q[i] - q_u[i]) * (q[i] - q_u[i]);
        norm_vec_q_u += (q_u[i]) * (q_u[i]);
    }
    cout << endl;
    cout << "Относительная норма вектора погрешности полученного решения:" << endl;
    cout << sqrt(norm_vec_err) / sqrt(norm_vec_q_u) << endl;
}
#pragma endregion


int main()
{
    setlocale(LC_ALL, "Russian");
    Create_Grid();
    BuildMatrA();
    BuildVecB();
    ConsiderBoundConditSecType();
    ConsiderBoundConditFirstType();
    ConsiderFictitiousNodes();
    m = NUM__NODES_IN_EDGE_X - 2;
    q.resize(NUM_NODES, 0);
    int max_iter = 1000;
    double eps = 1e-15;
    double w = 1;
    GaussZaid(w, eps, max_iter);
    Test();
    Output();
}
