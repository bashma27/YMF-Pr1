#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <set>
#include <iomanip>
#include <math.h>

using namespace std;
#pragma region Объявление глобальных переменных
int num_split_edge, num_nodes, num_split, choice; //число узлов в ребре, число узлов, количество разбиений, выбор типа сетки
double h_x, h_y; // шаг разбиений по x и y
double coef_q; // коэф. разрядки q
double lambda = 1; // коеф. лямбда
double gamma = 1; // коеф. гамма
int m; // расстояние между диагоналями
vector<vector<double>> al;
vector<double> di;
//vector<int> ia, ja, choice;
vector<double> q;
//vector<function<double(double, double)>> basic_func, deriv_basic_func_xi, deriv_basic_func_eta;
vector<pair<double, double>> nodes; // узлы сетки, заданные координатами 
vector<vector<vector<int>>> edge; // массив ребер с краевыми соответствующего рода

#pragma endregion

#pragma region Функции краевых условий и задание типа сетки
void EnterGridParameters() { // ввод параметров сетки
    cout << "Введите тип сетки:" << endl;
    cout << "1 - Равномерная;" << endl;
    cout << "2 - Неравномерная;" << endl << endl;
    cin >> choice;
    cout << "Введите количество разбиений ребра:" << endl;
    cin >> num_split;
    if (choice == 2) {
        cout << "Введите коэф. разрядки:" << endl;
        cin >> coef_q;
    }
    cout << endl;
}

double u_g(int var, double r, double z) { // краевое условие первого рода
    if (var == 0) {
        return exp(r * z);
    }
    else return exp(r * z);
}
double theta(int var, double r, double z) { // краевое условие второго рода
    if (var == 0) {
        return exp(r * z) * z;
    }
    else {
        return exp(r * z) * z;
    }
}

//double u_beta(int var, double r, double z) { // краевое условие третьего рода
//    if (var == 0) {
//        for (int i = 0; i < 4; i++) {
//            if (choice[i] != 3) continue;
//            return (exp(r * z) * z * normal[i * 2] + exp(r * z) * r * normal[i * 2 + 1]) + exp(r * z);
//        }
//    }
//    else {
//        bool first_condit = 1;
//        for (int i = 0; i < 4; i++) {
//            if (choice[i] != 3) continue;
//            if (first_condit) {
//                first_condit = 0;
//                continue;
//            }
//            return (exp(r * z) * z * normal[i * 2] + exp(r * z) * r * normal[i * 2 + 1]) + exp(r * z);
//        }
//    }
//}
#pragma endregion

#pragma region Задание функции
double f(double r, double z) {
    return exp(r * z) * (-z / r - z * z - r * r + 1 / r * r);
}

#pragma endregion

#pragma region Генерация сетки
void GenEndElGrid() { // создание конечноэлементной сетки
    ifstream file_in("grid_coordinates.txt");
    num_split_edge = num_split + 1; //число узлов в ребре 
    num_nodes = pow(num_split_edge, 2); //число узлов    
    nodes.resize(num_nodes);
    m = num_split;
    double start_point_x, start_point_y, end_point_x, end_point_y;
    file_in >> start_point_x >> end_point_x;
    file_in >> start_point_y >> end_point_y;
    file_in.close();
    if (choice == 1) {
        h_x = (end_point_x - start_point_x) / num_split;
        h_y = (end_point_y - start_point_y) / num_split;
        for (int i = 0, k = 0; i < num_nodes;) {
            for (int j = 0; j < num_split_edge; j++) {
                nodes[i + j] = { start_point_x + j * h_x, start_point_y + k * h_y };
            }
            i += num_split_edge;
            k++;
        }
    }
    else {   
        for (int j = 0; j < num_split_edge; j++) {
            if (j == 0) {
                h_x = 0;
                nodes[j] = { start_point_x + h_x, start_point_y };
                h_x = (end_point_x - start_point_x) * (1 - coef_q) / (1 - pow(coef_q, num_split));
            }
            else {
                nodes[j] = { nodes[j - 1].first + h_x, start_point_y};
                h_x *= coef_q;
            }
        }
        h_y = (end_point_y - start_point_y) * (1 - coef_q) / (1 - pow(coef_q, num_split));
        for (int i = num_split_edge; i < num_nodes;) {
            for (int j = 0; j < num_split_edge; j++) {
                if (j == 0) {
                    h_x = 0;
                    nodes[i + j] = { start_point_x + h_x, nodes[i - 1].second + h_y };
                    h_x = (end_point_x - start_point_x) * (1 - coef_q) / (1 - pow(coef_q, num_split));
                }
                else {
                    nodes[i + j] = { nodes[j - 1].first + h_x, nodes[i - 1].second + h_y };
                    h_x *= coef_q;
                }
            }
            h_y *= coef_q;
            i += num_split_edge;
        }
    }
}
#pragma endregion

#pragma region Работа с краевыми условиями

void BoundCondit() { // краевые условия
    edge.resize(2);
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
    }  
}

void ConsiderBoundConditFirstType(int node_num) { // учет краевых условий первого типа
   


}

void ConsiderBoundConditSecType() { // учет краевых условий второго типа
    for (int i = 1; i < num_split_edge - 1; i++) {
        int up_point = i + num_split_edge;
        if (choice == 1) {
            al[0][up_point] = -lambda / h_y;
            di[i] = lambda / h_y;
        }
        else {
            al[0][up_point] = -lambda / (nodes[up_point].second - nodes[i].second);
            di[i] = lambda / (nodes[up_point].second - nodes[i].second);
        }
    }
}
#pragma endregion

#pragma region Построение глобальной матрицы и вектора

bool FindInd(int i) {
    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < num_split_edge; j++) {
            if (edge[0][k][j] == i) {
                return true;
            }
        }
    }
    for (int j = 0; j < num_split_edge; j++) {
        if (edge[1][0][j] == i) {
            return true;
        }
    }
    return false;
}

void BuildMatrA() {
    al.assign(2, vector<double> (num_nodes));
    di.resize(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        if (FindInd(i)) continue;

        if (choice == 1) {     
            al[0][i] = - lambda / pow(h_y, 2);
            al[1][i] = - lambda / pow(h_x, 2);
            di[i] = 2 * lambda * (1 / pow(h_x, 2) + 1 / pow(h_y, 2)) + gamma;
        }

        else {
            int left_point, right_point, down_point, up_point;
            double h_x_prev, h_x_curr, h_y_prev, h_y_curr; // предыдущий / текущий
            down_point = i - num_split_edge;
            up_point = i + num_split_edge;
            left_point = i - 1;
            right_point = i + 1;
            h_x_prev = nodes[i].first - nodes[left_point].first;
            h_x_curr = nodes[right_point].first - nodes[i].first;
            h_y_prev = nodes[i].second - nodes[down_point].second;
            h_y_curr = nodes[up_point].second - nodes[i].second;
            al[0][i] = - 2 * lambda / (h_y_prev * (h_y_curr + h_y_prev));
            al[1][i] = - 2 * lambda / (h_x_prev * (h_x_curr + h_x_prev));
            di[i] = 2 * lambda * (1 / (h_x_curr * h_x_prev) + 1 / (h_y_curr * h_y_prev)) + gamma;
        }

    }
    ConsiderBoundConditSecType();
}

void BuildVecB() {
    //for (int k = 0; k < num_split * num_split; k++) {
    //    vector<int> node_num(9);
    //    for (int i = 0; i < 9; i++) {
    //        node_num[i] = array_quad[k][i];
    //    }
    //    vector<double> _r(9), _z(9);
    //    for (int i = 0; i < 9; i++) {
    //        _r[i] = nodes[node_num[i]].first, _z[i] = nodes[node_num[i]].second;
    //    }
    //    vector<double> a(3);
    //    a[0] = (_r[2] - _r[0]) * (_z[6] - _z[0]) - (_z[2] - _z[0]) * (_r[6] - _r[0]);
    //    a[1] = (_r[2] - _r[0]) * (_z[8] - _z[6]) - (_z[2] - _z[0]) * (_r[8] - _r[6]);
    //    a[2] = (_r[8] - _r[2]) * (_z[6] - _z[0]) - (_z[8] - _z[2]) * (_r[6] - _r[0]);
    //    vector<double> p(9);
    //    for (int i = 0; i < 9; i++) {
    //        p[i] = f(_r[i], _z[i]);
    //    }

    //    int sign_a0;
    //    if (a[0] > 0) sign_a0 = 1;
    //    else sign_a0 = -1;

    //    localVecB(_r, _z, p, a, node_num, sign_a0);
    //}
}
#pragma endregion

#pragma region Тестирование
void Test() {
    vector<double> q_u(num_nodes, 0);
    for (int i = 0; i < num_nodes; i++) {
        q_u[i] = u_g(0, nodes[i].first, nodes[i].second);
    }
    double norm_vec_err = 0, norm_vec_q_u = 0; // норма вектора погрешности и q_u
    for (int i = 0; i < num_nodes; i++) {
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
    EnterGridParameters();
    GenEndElGrid();
    BoundCondit();
    BuildMatrA();
    /*
    BuildVecB();
    ConsiderBoundCondit();
    q.resize(num_nodes, 0);
    vector<double> r(num_nodes);
    vector<double> z(num_nodes);
    vector<double> Mult(num_nodes);
    vector<double> Az(num_nodes);
    int max_iter = 1000;
    double eps = 1e-15;
    Test();*/
}
