#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <set>
#include <iomanip>
#include <math.h>
#define EPS 1e-16 //Для сравнения double

using namespace std;

#pragma region Объявление глобальных переменных
int num_split_edge_x, num_split_edge_y, num_nodes; //число узлов в ребре, число узлов
int Nx, Ny, L; // Nx - число вертикальных границ, Ny - число горизонтальных границ, L - число прямоугольников
//double h_x, h_y; // шаг разбиений по x и y
double lambda = 1; // коеф. лямбда
double gamma = 1; // коеф. гамма
int num_bc1, num_bc2; //Количество рёбер с первыми и вторыми краевыми условиями
int m; // расстояние между диагоналями
vector<vector<double>> al, au;
vector<double> Xw, Yw; // границы области
vector<vector<int>> W; // прямоугольники
vector<double> di, q, b;
vector<pair<double, double>> nodes; // узлы сетки, заданные координатами 
vector<vector<int>> edge; // первый индекс номер краевого, потом индексы
vector<vector<pair<int, double>>> koef; // Коэффициеты n, q  для каждой области в виде пары

#pragma endregion

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
            return true;
    }
    return false; // если по результату прохода узел не находится не в одном прямоугольнике то фиктивный узел.
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
    num_split_edge_x = CalcSum_ni(0, 0, Nx - 1) + 1;
    num_split_edge_y = CalcSum_ni(1, 0, Ny - 1) + 1;
    num_nodes = num_split_edge_x * num_split_edge_y;
    nodes.resize(num_nodes);
    return 0;
}

bool InVector(int value, vector<int>& temp) { //Для записи в более коротком формате.
    for (auto var : temp) {
        if (var == value)
            return false;
    }
    return true;
}

#pragma region Задание функций

double f(double x, double y) {
    return y * x;
    //return y * y - 2;
}

double theta(double x, double y) { // краевое условие второго рода
    return - x;
    //return - 2 * y;
}

double u_g(double x, double y) { // краевое условие первого рода
    return y * x;
    //return y * y;
}

#pragma endregion

#pragma region Вывод вектора решений
void Output() {
    ofstream f_result;
    f_result.imbue(locale("Russian"));
    f_result.open("f_result.txt");
    for (int i = 0; i < num_nodes; i++) {
        f_result << setprecision(15) << q[i] << endl;
    }
    f_result.close();
}
#pragma endregion

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
    File_Edge >> num_bc1 >> num_bc2; //считывание кол-во ребер каждого краевого условия
    edge.resize(2);
    vector<int> corner(num_bc1 + num_bc2, -1);
    it = 0;
    for (int l = 0; l < num_bc1; l++) { // Сначало обрабатываем 1 краевые
        File_Edge >> x1 >> x2 >> y1 >> y2;
        // Находим номер первого узла ребра (для наглядности разделил на две операции по каждому индексу)
        i = CalcSum_ni(0, 0, x1 - 1);
        j = CalcSum_ni(1, 0, y1 - 1) * num_split_edge_x;
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
                edge[0].push_back(ij + k * num_split_edge_x);
            if (InVector(ij + ny * num_split_edge_x, corner)) {
                edge[0].push_back(ij + ny * num_split_edge_x);
                corner[it] = ij + ny * num_split_edge_x;
                it++;
            }
        }
    }
    for (int l = 0; l < num_bc2; l++) { //Обрабатываем вторые краевыеы
        File_Edge >> x1 >> x2 >> y1 >> y2;
        //Находим номер первого узла ребра (для наглядности разделил на две операции по каждому индексу)
        i = CalcSum_ni(0, 0, x1 - 1);
        j = CalcSum_ni(1, 0, y1 - 1) * num_split_edge_x;
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
                edge[1].push_back(ij + k * num_split_edge_x);
            if (InVector(ij + ny * num_split_edge_x, corner)) {
                edge[1].push_back(ij + ny * num_split_edge_x);
                corner[it] = ij + ny * num_split_edge_x;
                it++;
            }
        }
    }
    File_Edge.close();
    corner.clear();
}

#pragma region Генерация сетки
//void GenEndElGrid() { // создание конечноэлементной сетки
//    ifstream file_in("grid_coordinates.txt");
//    num_split_edge = num_split + 1; //число узлов в ребре 
//    num_nodes = pow(num_split_edge, 2); //число узлов    
//    nodes.resize(num_nodes);
//    m = num_split - 1;
//    double start_point_x, start_point_y, end_point_x, end_point_y;
//    file_in >> start_point_x >> end_point_x;
//    file_in >> start_point_y >> end_point_y;
//    file_in.close();
//    if (choice == 1) {
//        h_x = (end_point_x - start_point_x) / num_split;
//        h_y = (end_point_y - start_point_y) / num_split;
//        for (int i = 0, k = 0; i < num_nodes;) {
//            for (int j = 0; j < num_split_edge; j++) {
//                nodes[i + j] = { start_point_x + j * h_x, start_point_y + k * h_y };
//            }
//            i += num_split_edge;
//            k++;
//        }
//    }
//    else {   
//        for (int j = 0; j < num_split_edge; j++) {
//            if (j == 0) {
//                h_x = 0;
//                nodes[j] = { start_point_x + h_x, start_point_y };
//                h_x = (end_point_x - start_point_x) * (1 - coef_q_x) / (1 - pow(coef_q_x, num_split));
//            }
//            else {
//                nodes[j] = { nodes[j - 1].first + h_x, start_point_y};
//                h_x *= coef_q_x;
//            }
//        }
//        h_y = (end_point_y - start_point_y) * (1 - coef_q_y) / (1 - pow(coef_q_y, num_split));
//        for (int i = num_split_edge; i < num_nodes;) {
//            for (int j = 0; j < num_split_edge; j++) {
//                if (j == 0) {
//                    h_x = 0;
//                    nodes[i + j] = { start_point_x + h_x, nodes[i - 1].second + h_y };
//                    h_x = (end_point_x - start_point_x) * (1 - coef_q_x) / (1 - pow(coef_q_x, num_split));
//                }
//                else {
//                    nodes[i + j] = { nodes[j - 1].first + h_x, nodes[i - 1].second + h_y };
//                    h_x *= coef_q_x;
//                }
//            }
//            h_y *= coef_q_y;
//            i += num_split_edge;
//        }
//    }
//}
#pragma endregion

#pragma region Работа с краевыми условиями

void BoundCondit() { // краевые условия
/*    edge.resize(2);
    edge[0].resize(3);
    int start_point = num_split_edge - 1;
    for (int i = 0; i < num_nodes;) {
        edge[0][0].push_back(start_point + i);
        i += num_split_edge;
    }
    start_point = num_nodes - num_split_edge;
    for (int i = 0; i < num_split_edge; i++) {
        edge[0][1].push_back(start_point + i);
    }
    start_point = 0;
    for (int i = 0; i < num_nodes;) {
        edge[0][2].push_back(start_point + i);
        i += num_split_edge;
    }
    edge[1].resize(1);
    for (int i = 0; i < num_split_edge; i++) {
        edge[1][0].push_back(start_point + i);
    } */ 
}

//void ConsiderBoundConditFirstType() { // учет краевых условий первого типа
//   
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < num_split_edge; j++) {
//            int _i = edge[0][i][j];
//            b[_i] = u_g(nodes[_i].first, nodes[_i].second);
//            di[_i] = double(1);           
//        }
//    }
//
//}
//
//void ConsiderBoundConditSecType() { // учет краевых условий второго типа
//    for (int i = 1; i < num_split_edge - 1; i++) {
//        int up_point = i + num_split_edge;
//        if (choice == 1) {
//            au[0][up_point] = -lambda / h_y;
//            di[i] = lambda / h_y;
//        }
//        else {
//            au[0][up_point] = -lambda / (nodes[up_point].second - nodes[i].second);
//            di[i] = lambda / (nodes[up_point].second - nodes[i].second);
//        }
//        b[i] = theta(nodes[i].first, nodes[i].second);
//    }
//}
//#pragma endregion
//
//#pragma region Построение глобальной матрицы и вектора
//
//bool FindInd(int i) {
//    for (int k = 0; k < 3; k++) {
//        for (int j = 0; j < num_split_edge; j++) {
//            if (edge[0][k][j] == i) {
//                return true;
//            }
//        }
//    }
//    for (int j = 0; j < num_split_edge; j++) {
//        if (edge[1][0][j] == i) {
//            return true;
//        }
//    }
//    return false;
//}
//
//void BuildMatrA() {
//    al.assign(2, vector<double> (num_nodes));
//    au.assign(2, vector<double> (num_nodes));
//    di.resize(num_nodes);
//    for (int i = 0; i < num_nodes; i++) {
//        if (FindInd(i)) continue;
//
//        if (choice == 1) {     
//            al[0][i] = - lambda / pow(h_y, 2);
//            al[1][i] = - lambda / pow(h_x, 2);
//            au[0][i + num_split + 1] = - lambda / pow(h_y, 2);
//            au[1][i + 1] = - lambda / pow(h_x, 2);
//            di[i] = 2 * lambda * (1 / pow(h_x, 2) + 1 / pow(h_y, 2)) + gamma;
//        }
//
//        else {
//            int left_point, right_point, down_point, up_point;
//            double h_x_prev, h_x_curr, h_y_prev, h_y_curr; // предыдущий / текущий
//            down_point = i - num_split_edge;
//            up_point = i + num_split_edge;
//            left_point = i - 1;
//            right_point = i + 1;
//            h_x_prev = nodes[i].first - nodes[left_point].first;
//            h_x_curr = nodes[right_point].first - nodes[i].first;
//            h_y_prev = nodes[i].second - nodes[down_point].second;
//            h_y_curr = nodes[up_point].second - nodes[i].second;
//            al[0][i] = - 2 * lambda / (h_y_prev * (h_y_curr + h_y_prev));
//            al[1][i] = - 2 * lambda / (h_x_prev * (h_x_curr + h_x_prev));
//            au[0][i + num_split + 1] = - 2 * lambda / (h_y_curr * (h_y_curr + h_y_prev));
//            au[1][i + 1] = - 2 * lambda / (h_x_curr * (h_x_curr + h_x_prev));
//            di[i] = 2 * lambda * (1 / (h_x_curr * h_x_prev) + 1 / (h_y_curr * h_y_prev)) + gamma;
//        }
//    }
//}

//void BuildVecB() {
//    b.resize(num_nodes);
//    for (int i = 0; i < num_nodes; i++) {
//        if (FindInd(i)) continue;
//        b[i] = f(nodes[i].first, nodes[i].second);
//    }
//}

#pragma endregion

#pragma region Решение СЛАУ(Гаусс-Зейдель)

double norma() {
    double res = 0;
    for (int i = 0; i < num_nodes; i++) {
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

    if (i < num_nodes - 1) {
        temp += au[1][i + 1] * q[i + 1];
    }
    if (i < num_nodes - m - 2) {
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
        for (int i = 0; i < num_nodes; i++) {
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
//void Test() {
//    vector<double> q_u(num_nodes, 0);
//    for (int i = 0; i < num_nodes; i++) {
//        q_u[i] = u_g(0, nodes[i].first, nodes[i].second);
//    }
//    double norm_vec_err = 0, norm_vec_q_u = 0; // норма вектора погрешности и q_u
//    for (int i = 0; i < num_nodes; i++) {
//        norm_vec_err += (q[i] - q_u[i]) * (q[i] - q_u[i]);
//        norm_vec_q_u += (q_u[i]) * (q_u[i]);
//    }
//    cout << endl;
//    cout << "Относительная норма вектора погрешности полученного решения:" << endl;
//    cout << sqrt(norm_vec_err) / sqrt(norm_vec_q_u) << endl;
//}
#pragma endregion

int main()
{
    setlocale(LC_ALL, "Russian");
    Create_Grid();
    //GenEndElGrid();
    //BoundCondit();
    //BuildMatrA();
    //BuildVecB();
    //ConsiderBoundConditSecType();
    //ConsiderBoundConditFirstType();   
    q.resize(num_nodes, 0);
    int max_iter = 1000;
    double eps = 1e-15;
    double w = 1;
    //GaussZaid(w, eps, max_iter);
    //Output();
}
