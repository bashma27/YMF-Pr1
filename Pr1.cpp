#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <set>
#include <iomanip>
#include <math.h>

using namespace std;
#pragma region Объявление глобальных переменных
int num_split_edge, num_nodes, num_split; //число узлов в ребре, число узлов, количество разбиений  
//double beta = 1, coef_z = 1.15, coef_r = 1.15; // параметр бета(константа), коэф. разрядки z, коэф. разрядки r
//vector<int> ia, ja, choice;
//vector<double> aal, di, b, q, L_sq, di_sq, normal;
//vector<function<double(double, double)>> basic_func, deriv_basic_func_xi, deriv_basic_func_eta;
vector<pair<double, double>> nodes; // узлы сетки, заданные координатами 
vector<vector<vector<int>>> edge; // массив ребер с краевыми соответствующего рода
vector<vector<int>> array_quad; // массив прямоугольников
//vector<vector<double>> A_1, A_2, A_3;

#pragma endregion

#pragma region Функции краевых условий и задание этих функций
//void EnterBoundCondit() { // ввод краевых условий
//    cout << "Введите тип краевого условия на соответствующем ребре:" << endl;
//    cout << "1 - первое краевое условие;" << endl;
//    cout << "2 - второе краевое условие;" << endl;
//    cout << "3 - третье краевое условие;" << endl << endl;
//    choice.resize(4);
//    cout << "Краевое условие нижнего ребра:" << endl;
//    cin >> choice[0];
//    cout << endl;
//    cout << "Краевое условие правого ребра:" << endl;
//    cin >> choice[1];
//    cout << endl;
//    cout << "Краевое условие верхнего ребра:" << endl;
//    cin >> choice[2];
//    cout << endl;
//    cout << "Краевое условие левого ребра:" << endl;
//    cin >> choice[3];
//    cout << endl;
//}

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

#pragma region Задание функции, лямбды и якобиана преобразования
double f(double r, double z) {
    return exp(r * z) * (-z / r - z * z - r * r + 1 / r * r);
}

double lambda(double z) {
    return 1;
}
#pragma endregion

#pragma region Генерация сетки
void GenEndElGrid() { // создание конечноэлементной сетки
    ifstream file_in("grid_coordinates.txt");
    file_in >> num_split;
    num_split_edge = 2 * num_split + 1; //число узлов в ребре 
    num_nodes = pow(num_split_edge, 2); //число узлов    
    nodes.resize(num_nodes);
    vector<pair<double, double>> step(4); // шаг разбиений на конкретном ребре
    double a, b;
    int l_d, r_d, l_u, r_u; //левый нижный, правный нижний, правый верхний, левый верхний соответственно   
    l_d = 0;
    r_d = num_split_edge - 1;
    l_u = num_nodes - num_split_edge;
    r_u = num_nodes - 1;
    file_in >> a >> b;
    nodes[l_d] = { a, b };
    file_in >> a >> b;
    nodes[r_d] = { a, b };
    file_in >> a >> b;
    nodes[r_u] = { a, b };
    file_in >> a >> b;
    nodes[l_u] = { a, b };
    file_in.close();
    int temp_1, temp_2;
    if (coef_r > 1) temp_1 = 1 / (coef_r - 1);
    else temp_1 = 1 / coef_r;
    if (coef_z > 1) temp_2 = 1 / (coef_z - 1);
    else temp_2 = 1 / coef_z;
    int num_eq_parts_r, num_eq_parts_z;
    if (coef_r != 1) {
        num_eq_parts_r = (1 - pow(temp_1, num_split)) / (1 - temp_1); // количество равных частей(по формуле геом. прог.)
    }
    else num_eq_parts_r = num_split;
    if (coef_z != 1) {
        num_eq_parts_z = (1 - pow(temp_2, num_split)) / (1 - temp_2); // количество равных частей(по формуле геом. прог.)
    }
    else num_eq_parts_z = num_split;

    step[0] = { (nodes[r_d].first - nodes[l_d].first) / num_eq_parts_r, (nodes[r_d].second - nodes[l_d].second) / num_eq_parts_z }; //шаг(разбиение) для нижнего ребра
    step[1] = { (nodes[r_u].first - nodes[r_d].first) / num_eq_parts_r, (nodes[r_u].second - nodes[r_d].second) / num_eq_parts_z }; //шаг(разбиение) для правого ребра
    step[2] = { (nodes[r_u].first - nodes[l_u].first) / num_eq_parts_r, (nodes[r_u].second - nodes[l_u].second) / num_eq_parts_z }; //шаг(разбиение) для верхнего ребра
    step[3] = { (nodes[l_u].first - nodes[l_d].first) / num_eq_parts_r, (nodes[l_u].second - nodes[l_d].second) / num_eq_parts_z }; //шаг(разбиение) для левого ребра
    if (coef_r > 1) {
        for (int i = 0; i < 4; i++) {
            step[i].first *= temp_1;
        }
    }
    if (coef_z > 1) {
        for (int i = 0; i < 4; i++) {
            step[i].second *= temp_2;
        }
    }
    for (int j = 2; j < num_split_edge - 1; j++) {
        if (j % 2 != 0) continue;
        nodes[j].first = nodes[j - 2].first + step[0].first; //нефиктивные узлы нижнего ребра
        nodes[j].second = (nodes[r_d].second - nodes[l_d].second) * ((nodes[j].first - nodes[l_d].first) / (nodes[r_d].first - nodes[l_d].first)) + nodes[l_d].second; //нефиктивные узлы нижнего ребра
        if (coef_r > 1) step[0].first /= temp_1;
        else step[0].first *= temp_1;
    }
    for (int j = 1; j < num_split_edge - 1; j++) {
        if (j % 2 == 0) continue;
        nodes[j] = { (nodes[j - 1].first + nodes[j + 1].first) / 2, (nodes[j - 1].second + nodes[j + 1].second) / 2 }; //фиктивные узлы нижнего ребра
    }
    for (int j = 2; j < num_split_edge - 1; j++) {
        if (j % 2 != 0) continue;
        nodes[l_u + j].first = nodes[l_u + j - 2].first + step[2].first; //нефиктивные узлы верхнего ребра
        nodes[l_u + j].second = (nodes[r_u].second - nodes[l_u].second) * ((nodes[l_u + j].first - nodes[l_u].first) / (nodes[r_u].first - nodes[l_u].first)) + nodes[l_u].second; //нефиктивные узлы верхнего ребра
        if (coef_r > 1) step[2].first /= temp_1;
        else step[2].first *= temp_1;
    }
    for (int j = 1; j < num_split_edge - 1; j++) {
        if (j % 2 == 0) continue;
        nodes[l_u + j] = { (nodes[l_u + j - 1].first + nodes[l_u + j + 1].first) / 2, (nodes[l_u + j - 1].second + nodes[l_u + j + 1].second) / 2 }; //фиктивные узлы верхнего ребра
    }
    for (int i = 0; i < num_nodes;) {
        if (i == l_d || i == l_u) {
            i += num_split_edge;
            continue;
        }
        if (i % 2 * num_split_edge != 0) {
            i += num_split_edge;
            continue;
        }
        double in_step_r;
        nodes[i].second = nodes[i - 2 * num_split_edge].second + step[3].second; //начальная вершина в строке по z
        nodes[i].first = (nodes[l_u].first - nodes[l_d].first) * ((nodes[i].second - nodes[l_d].second) / (nodes[l_u].second - nodes[l_d].second)) + nodes[l_d].first; //начальная вершина в строке по r
        nodes[i + num_split_edge - 1].second = nodes[i - num_split_edge - 1].second + step[1].second; //конечная вершина в строке по z
        nodes[i + num_split_edge - 1].first = (nodes[r_u].first - nodes[r_d].first) * ((nodes[i + num_split_edge - 1].second - nodes[r_d].second) / (nodes[r_u].second - nodes[r_d].second)) + nodes[r_d].first; //конечная вершина в строке по r
        in_step_r = (nodes[i + num_split_edge - 1].first - nodes[i].first) / num_eq_parts_r;
        if (coef_r > 1) in_step_r *= temp_1;
        for (int j = 2; j < num_split_edge - 1; j++) {
            if (j % 2 != 0) continue;
            nodes[i + j].first = nodes[i + j - 2].first + in_step_r;
            nodes[i + j].second = (nodes[i + num_split_edge - 1].second - nodes[i].second) * ((nodes[i + j].first - nodes[i].first) / (nodes[i + num_split_edge - 1].first - nodes[i].first)) + nodes[i].second;
            if (coef_r > 1) in_step_r /= temp_1;
            else in_step_r *= temp_1;
        }
        for (int j = 1; j < num_split_edge - 1; j++) {
            if (j % 2 == 0) continue;
            nodes[i + j] = { (nodes[i + j - 1].first + nodes[i + j + 1].first) / 2, (nodes[i + j - 1].second + nodes[i + j + 1].second) / 2 };
        }
        i += num_split_edge;
        if (coef_z > 1) step[3].second /= temp_2, step[1].second /= temp_2;
        else step[3].second *= temp_2, step[1].second *= temp_2;
    }
    for (int i = 0; i < num_nodes;) {
        if (i == l_d || i == l_u) {
            i += num_split_edge;
            continue;
        }
        if (i % 2 * num_split_edge == 0) {
            i += num_split_edge;
            continue;
        }
        double in_step_r;
        nodes[i] = { (nodes[i + num_split_edge].first + nodes[i - num_split_edge].first) / 2, (nodes[i + num_split_edge].second + nodes[i - num_split_edge].second) / 2 }; //начальная вершина в строке
        nodes[i + num_split_edge - 1] = { (nodes[i + 2 * num_split_edge - 1].first + nodes[i - 1].first) / 2, (nodes[i + 2 * num_split_edge - 1].second + nodes[i - 1].second) / 2 }; //конечная вершина в строке
        in_step_r = (nodes[i + num_split_edge - 1].first - nodes[i].first) / num_eq_parts_r;
        if (coef_r > 1) in_step_r *= temp_1;
        for (int j = 2; j < num_split_edge - 1; j++) {
            if (j % 2 != 0) continue;
            nodes[i + j].first = nodes[i + j - 2].first + in_step_r;
            nodes[i + j].second = (nodes[i + num_split_edge - 1].second - nodes[i].second) * ((nodes[i + j].first - nodes[i].first) / (nodes[i + num_split_edge - 1].first - nodes[i].first)) + nodes[i].second;
            if (coef_r > 1) in_step_r /= temp_1;
            else in_step_r *= temp_1;
        }
        for (int j = 1; j < num_split_edge - 1; j++) {
            if (j % 2 == 0) continue;
            nodes[i + j] = { (nodes[i + j - 1].first + nodes[i + j + 1].first) / 2, (nodes[i + j - 1].second + nodes[i + j + 1].second) / 2 };
        }
        i += num_split_edge;
    }
}
#pragma endregion

#pragma region Создание массива конечных элементов 
//void ArrayQuadrilaterals() { // массив четырехугольников
//    array_quad.resize(num_split * num_split);
//    int k = 0;
//    for (int i = 0; i < num_split * num_split;) {
//        for (int j = 0; j < num_split; j++) {
//            int l_d_node = k + 2 * j; // нижний левый номер узла конечного элемента   
//            array_quad[i + j].push_back(l_d_node);
//            array_quad[i + j].push_back(l_d_node + 1);
//            array_quad[i + j].push_back(l_d_node + 2);
//            array_quad[i + j].push_back(l_d_node + num_split_edge);
//            array_quad[i + j].push_back(l_d_node + 1 + num_split_edge);
//            array_quad[i + j].push_back(l_d_node + 2 + num_split_edge);
//            array_quad[i + j].push_back(l_d_node + 2 * num_split_edge);
//            array_quad[i + j].push_back(l_d_node + 1 + 2 * num_split_edge);
//            array_quad[i + j].push_back(l_d_node + 2 + 2 * num_split_edge);
//        }
//        i += num_split;
//        k += 2 * num_split_edge;
//    }
//}
#pragma endregion

#pragma region Работа с краевыми условиями

int GetStartPoint(int i) { // задаем начальную точку ребра
    switch (i)
    {
    case 0:
        return 0;
        break;
    case 1:
        return num_split_edge - 1;
        break;
    case 2:
        return num_nodes - num_split_edge;
        break;
    case 3:
        return 0;
        break;
    }
}

int GetStepEndEl(int i) { // задаем шаг до след конеч. эл.
    if (i == 0 || i == 2) {
        return 2;
    }
    else return 2 * num_split_edge;
}

void BoundCondit() { // краевые условия
    edge.resize(3);
    int num_func_1 = 0, num_func_2 = 0, num_func_3 = 0;
    for (int i = 0; i < 4; i++) {
        int start_point = GetStartPoint(i);
        int step = GetStepEndEl(i);
        switch (choice[i])
        {
        case 1:
            for (int j = 0; j < num_split_edge; j++) {
                edge[0].resize(num_split_edge * (num_func_1 + 1));
                edge[0][j + num_split_edge * num_func_1].push_back(start_point + j * step / 2);
                edge[0][j + num_split_edge * num_func_1].push_back(num_func_1);
            }
            num_func_1++;
            break;
        case 2:
            for (int j = 0; j < num_split; j++) {
                edge[1].resize(num_split * (num_func_2 + 1));
                edge[1][j + num_split * num_func_2].push_back(start_point);
                edge[1][j + num_split * num_func_2].push_back(start_point + step / 2);
                edge[1][j + num_split * num_func_2].push_back(start_point + step);
                edge[1][j + num_split * num_func_2].push_back(num_func_2);
                start_point += step;
            }
            num_func_2++;
            break;
        case 3:
            for (int j = 0; j < num_split; j++) {
                edge[2].resize(num_split * (num_func_3 + 1));
                edge[2][j + num_split * num_func_3].push_back(start_point);
                edge[2][j + num_split * num_func_3].push_back(start_point + step / 2);
                edge[2][j + num_split * num_func_3].push_back(start_point + step);
                edge[2][j + num_split * num_func_3].push_back(num_func_3);
                start_point += step;
            }
            num_func_3++;
            break;
        }
    }
}

bool FindInd(int i) {
    for (int j = 0; j < edge[0].size(); j++) {
        if (edge[0][j][0] == i) {
            return true;
        }
    }
    return false;
}

void ConsiderBoundConditFirstType(int node_num) { // учет краевых условий первого типа
    double _r, _z;
    _r = nodes[edge[0][node_num][0]].first, _z = nodes[edge[0][node_num][0]].second;
    b[edge[0][node_num][0]] = u_g(edge[0][node_num][1], _r, _z);
    di[edge[0][node_num][0]] = 1;
    for (int i = ia[edge[0][node_num][0]]; i < ia[edge[0][node_num][0] + 1]; i++) {
        int _i = ja[i];
        if (FindInd(_i)) {
            aal[i] = 0;
            continue;
        }
        b[_i] -= b[edge[0][node_num][0]] * aal[i];
        aal[i] = 0;
    }
    for (int i = edge[0][node_num][0]; i < num_nodes; i++) {
        int k = 0;
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            if (ja[j] == edge[0][node_num][0]) {
                if (FindInd(i)) {
                    aal[j] = 0;
                    continue;
                }
                b[i] -= b[edge[0][node_num][0]] * aal[j];
                aal[j] = 0;
            }
        }
    }
}

void ConsiderBoundConditSecType(int num_edge) { // учет краевых условий второго типа
    vector<double> _r(3), _z(3);
    _r[0] = nodes[edge[1][num_edge][0]].first, _z[0] = nodes[edge[1][num_edge][0]].second;
    _r[1] = nodes[edge[1][num_edge][1]].first, _z[1] = nodes[edge[1][num_edge][1]].second;
    _r[2] = nodes[edge[1][num_edge][2]].first, _z[2] = nodes[edge[1][num_edge][2]].second;
    double h = sqrt(pow(_r[0] - _r[2], 2) + pow(_z[0] - _z[2], 2));
    vector<double> b_s2(3, 0);
    vector<double> _theta(3);
    _theta[0] = theta(edge[1][num_edge][3], _r[0], _z[0]);
    _theta[1] = theta(edge[1][num_edge][3], _r[1], _z[1]);
    _theta[2] = theta(edge[1][num_edge][3], _r[2], _z[2]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            b_s2[i] += (A_1[i][j] * _r[0] + A_2[i][j] * _r[1] + A_3[i][j] * _r[2]) * _theta[j];
        }
        b_s2[i] *= h * beta / double(420);
    }
    AddLocalVecBound(1, num_edge, b_s2);
}

//void ConsiderBoundConditThirdType(int num_edge) { // учет краевых условий третьего типа
//    vector<double> _r(3), _z(3);
//    _r[0] = nodes[edge[2][num_edge][0]].first, _z[0] = nodes[edge[2][num_edge][0]].second;
//    _r[1] = nodes[edge[2][num_edge][1]].first, _z[1] = nodes[edge[2][num_edge][1]].second;
//    _r[2] = nodes[edge[2][num_edge][2]].first, _z[2] = nodes[edge[2][num_edge][2]].second;
//    double h = sqrt(pow(_r[0] - _r[2], 2) + pow(_z[0] - _z[2], 2));
//    vector<vector<double>> _A(3);
//    _A[0].resize(1);
//    _A[1].resize(2);
//    _A[2].resize(3);
//    vector<double> b_s3(3, 0);
//    vector<double> _u_beta(3);
//    _u_beta[0] = u_beta(edge[2][num_edge][3], _r[0], _z[0]);
//    _u_beta[1] = u_beta(edge[2][num_edge][3], _r[1], _z[1]);
//    _u_beta[2] = u_beta(edge[2][num_edge][3], _r[2], _z[2]);
//    for (int i = 0; i < 3; i++) {
//        for (int j = i; j < 3; j++) {
//            _A[j][i] = h * beta / double(420) * (A_1[j][i] * _r[0] + A_2[j][i] * _r[1] + A_3[j][i] * _r[2]);
//        }
//    }
//    AddLocalMartBound_3(num_edge, _A);
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < i; j++) {
//            b_s3[i] += _A[i][j] * _u_beta[j];
//        }
//        for (int j = i; j < 3; j++) {
//            b_s3[i] += _A[j][i] * _u_beta[j];
//        }
//    }
//    AddLocalVecBound(2, num_edge, b_s3);
//}

void ConsiderBoundCondit() { // учет всех краевых
    //for (int i = 0; i < edge[2].size(); i++) { // учет третьих краевых
    //    ConsiderBoundConditThirdType(i);
    //}
    for (int i = 0; i < edge[1].size(); i++) { // учет вторых краевых
        ConsiderBoundConditSecType(i);
    }
    for (int i = 0; i < edge[0].size(); i++) { // учет первых краевых(не в верхнем цикле, тк должен быть в самом конце)
        ConsiderBoundConditFirstType(i);
    }
}
#pragma endregion

#pragma region Построение глобальной матрицы и вектора
void BuildMatrA() {
    for (int i = 0; i < num_split * num_split; i++) {
        vector<int> node_num(9);
        for (int j = 0; j < 9; j++) {
            node_num[j] = array_quad[i][j];
        }
        vector<double> _r(4), _z(4);
        _r[0] = nodes[node_num[0]].first, _z[0] = nodes[node_num[0]].second;
        _r[1] = nodes[node_num[2]].first, _z[1] = nodes[node_num[2]].second;
        _r[2] = nodes[node_num[6]].first, _z[2] = nodes[node_num[6]].second;
        _r[3] = nodes[node_num[8]].first, _z[3] = nodes[node_num[8]].second;
        vector<double> a(3);
        a[0] = (_r[1] - _r[0]) * (_z[2] - _z[0]) - (_z[1] - _z[0]) * (_r[2] - _r[0]);
        a[1] = (_r[1] - _r[0]) * (_z[3] - _z[2]) - (_z[1] - _z[0]) * (_r[3] - _r[2]);
        a[2] = (_r[3] - _r[1]) * (_z[2] - _z[0]) - (_z[3] - _z[1]) * (_r[2] - _r[0]);

        int sign_a0;
        if (a[0] > 0) sign_a0 = 1;
        else sign_a0 = -1;

        LocalMatrMass(_r, _z, a, node_num, sign_a0);
        LocalMatrStiff(_r, _z, a, node_num, sign_a0);
    }
}

void BuildVecB() {
    for (int k = 0; k < num_split * num_split; k++) {
        vector<int> node_num(9);
        for (int i = 0; i < 9; i++) {
            node_num[i] = array_quad[k][i];
        }
        vector<double> _r(9), _z(9);
        for (int i = 0; i < 9; i++) {
            _r[i] = nodes[node_num[i]].first, _z[i] = nodes[node_num[i]].second;
        }
        vector<double> a(3);
        a[0] = (_r[2] - _r[0]) * (_z[6] - _z[0]) - (_z[2] - _z[0]) * (_r[6] - _r[0]);
        a[1] = (_r[2] - _r[0]) * (_z[8] - _z[6]) - (_z[2] - _z[0]) * (_r[8] - _r[6]);
        a[2] = (_r[8] - _r[2]) * (_z[6] - _z[0]) - (_z[8] - _z[2]) * (_r[6] - _r[0]);
        vector<double> p(9);
        for (int i = 0; i < 9; i++) {
            p[i] = f(_r[i], _z[i]);
        }

        int sign_a0;
        if (a[0] > 0) sign_a0 = 1;
        else sign_a0 = -1;

        localVecB(_r, _z, p, a, node_num, sign_a0);
    }
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
    EnterBoundCondit();
    GenEndElGrid();
    ArrayQuadrilaterals();
    BasicFunc();
    DerivBasicFuncXi();
    DerivBasicFuncEta();
    BoundCondit();
    GetVecNormal();
    BuildBoundMatrices();
    GeneratePortrait();
    BuildMatrA();
    BuildVecB();
    ConsiderBoundCondit();
    q.resize(num_nodes, 0);
    vector<double> r(num_nodes);
    vector<double> z(num_nodes);
    vector<double> Mult(num_nodes);
    vector<double> Az(num_nodes);
    int max_iter = 1000;
    double eps = 1e-15;
    LU_sq_MSG(q, r, z, Az, Mult, num_nodes, eps, max_iter);
    Test();
}
