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
vector<double> di, q, b;
//vector<int> ia, ja, choice;
//vector<function<double(double, double)>> basic_func, deriv_basic_func_xi, deriv_basic_func_eta;
vector<pair<double, double>> nodes; // узлы сетки, заданные координатами 
vector<vector<vector<int>>> edge; // массив ребер с краевыми соответствующего рода

#pragma endregion

#pragma region Задание типа сетки
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

#pragma region Задание функций

double f(double x, double y) {
    return y;
}

//double theta(double x, double y) { // краевое условие второго рода
//    return - 2 * y;
//}

double u_g(double x, double y) { // краевое условие первого рода
    return y;
}

#pragma endregion

#pragma region Генерация сетки
void GenEndElGrid() { // создание конечноэлементной сетки
    ifstream file_in("grid_coordinates.txt");
    num_split_edge = num_split + 1; //число узлов в ребре 
    num_nodes = pow(num_split_edge, 2); //число узлов    
    nodes.resize(num_nodes);
    m = num_split - 1;
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
    edge[0].resize(4);
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
    //edge[1].resize(1);
    for (int i = 0; i < num_split_edge; i++) {
        edge[0][3].push_back(start_point + i);
    }  
}

void ConsiderBoundConditFirstType() { // учет краевых условий первого типа
   
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < num_split_edge; j++) {
            int _i = edge[0][i][j];
            b[_i] = u_g(nodes[_i].first, nodes[_i].second);
            di[_i] = double(1);           
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < num_split_edge; j++) {
            int _i = edge[0][i][j];
            int _i_prev, _i_next, outer_di_prev, outer_di_next; // индексы ближней и внешней диагоналей
            _i_prev = _i - 1;
            _i_next = _i + 1;
            outer_di_prev = _i - num_split - 1;
            outer_di_next = _i + num_split + 1;
            if (al[0][_i] != 0) {
                b[outer_di_prev] -= b[_i] * al[0][_i];
                al[0][_i] = 0;
            }
            if (al[1][_i] != 0) {
                b[_i_prev] -= b[_i] * al[1][_i];
                al[1][_i] = 0;
            }
            if (_i_next >= num_nodes) continue;
            if (al[1][_i_next] != 0) {
                b[_i_next] -= b[_i] * al[1][_i_next];
                al[1][_i_next] = 0;
            }
            if (outer_di_next >= num_nodes) continue;
            if (al[0][outer_di_next] != 0) {
                b[outer_di_next] -= b[_i] * al[0][outer_di_next];
                al[0][outer_di_next] = 0;
            }

        }
    }

}

//void ConsiderBoundConditSecType() { // учет краевых условий второго типа
//    for (int i = 1; i < num_split_edge - 1; i++) {
//        int up_point = i + num_split_edge;
//        if (choice == 1) {
//            al[0][up_point] = -lambda / h_y;
//            di[i] = lambda / h_y;
//        }
//        else {
//            al[0][up_point] = -lambda / (nodes[up_point].second - nodes[i].second);
//            di[i] = lambda / (nodes[up_point].second - nodes[i].second);
//        }
//        b[i] = theta(nodes[i].first, nodes[i].second);
//    }
//}
#pragma endregion

#pragma region Построение глобальной матрицы и вектора

bool FindInd(int i) {
    for (int k = 0; k < 4; k++) {
        for (int j = 0; j < num_split_edge; j++) {
            if (edge[0][k][j] == i) {
                return true;
            }
        }
    }
    /*for (int j = 0; j < num_split_edge; j++) {
        if (edge[1][0][j] == i) {
            return true;
        }
    }*/
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
            al[0][i + num_split + 1] = -lambda / pow(h_y, 2);
            al[1][i + 1] = -lambda / pow(h_x, 2);
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
}

void BuildVecB() {
    b.resize(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        if (FindInd(i)) continue;
        b[i] = f(nodes[i].first, nodes[i].second);
    }
}

#pragma endregion

#pragma region Решение СЛАУ(Гаусс-Зейдель)

double norma(vector<double> v, int n) {
    double res = 0;
    for (int i = 0; i < n; i++) {
        res += v[i] * v[i];
    }
    return sqrt(res);
}

double iter(vector<double> x0, vector<double> x, int n, int m, int i) {
    double temp = 0;
    temp += di[i] * x0[i];

    if (i > 0) {
        temp += al[1][i] * x0[i - 1];
    }
    if (i > m + 1) {
        temp += al[0][i] * x0[i - m - 2];
    }

    if (i < n - 1) {
        temp += al[1][i] * x0[i + 1];
    }
    if (i < n - m - 2) {
        temp += al[0][i] * x0[i + m + 2];
    }
    temp = b[i] - temp;

    return temp;
}

void GaussZaid(vector<double> q, double w, double eps, int n, int m, int max_iter) {
    int counter = 1;
    double resid = 1;
    double norm = norma(b, n);
    for (; counter < max_iter && resid > eps; counter++) {
        resid = 0;
        for (int i = 0; i < n; i++) {
            double temp = iter(q, q, n, m, i);
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
    EnterGridParameters();
    GenEndElGrid();
    BoundCondit();
    BuildMatrA();
    BuildVecB();
    //ConsiderBoundConditSecType();
    ConsiderBoundConditFirstType();
    
    q.resize(num_nodes, 0);
    int max_iter = 1000;
    double eps = 1e-15;
    double w = 1;
    GaussZaid(q, w, eps, num_nodes, m, max_iter);
}
