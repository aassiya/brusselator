#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

const double A = 1.0;
const double B = 1.5;

// функция, описывающая систему Брюсселятора
void brusselator(double t, const std::vector<double> &z, std::vector<double> &dzdt) {
    double x = z[0];
    double y = z[1];

    dzdt[0] = A - (B + 1) * x + x * x * y; // dx/dt
    dzdt[1] = B * x - x * x * y;           // dy/dt
}

// метод Рунге-Кутты 4-го порядка
void rungeKutta4(double t, double h, std::vector<double> &z) {
    std::vector<double> k1(2), k2(2), k3(2), k4(2), z_temp(2);

    // вычисление k1
    brusselator(t, z, k1);

    // вычисление k2
    for (size_t i = 0; i < z.size(); ++i) {
        z_temp[i] = z[i] + 0.5 * h * k1[i];
    }
    brusselator(t + 0.5 * h, z_temp, k2);

    // вычисление k3
    for (size_t i = 0; i < z.size(); ++i) {
        z_temp[i] = z[i] + 0.5 * h * k2[i];
    }
    brusselator(t + 0.5 * h, z_temp, k3);

    // вычисление k4
    for (size_t i = 0; i < z.size(); ++i) {
        z_temp[i] = z[i] + h * k3[i];
    }
    brusselator(t + h, z_temp, k4);

    // обновление значений z
    for (size_t i = 0; i < z.size(); ++i) {
        z[i] += (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}

// метод Дорманда-Принса 5-го порядка
void dormandPrince5(double t, double &h, vector<double> &z, double tolerance) {
    std::vector<double> k1(2), k2(2), k3(2), k4(2), k5(2), k6(2), k7(2), z_temp(2);

    // коэффициенты метода
    const double a2 = 1.0 / 5.0;
    const double a3 = 3.0 / 10.0;
    const double a4 = 4.0 / 5.0;
    const double a5 = 8.0 / 9.0;
    const double a6 = 1.0;
    const double a7 = 1.0;

    const double b21 = 1.0 / 5.0;
    const double b31 = 3.0 / 40.0;
    const double b32 = 9.0 / 40.0;
    const double b41 = 44.0 / 45.0;
    const double b42 = -56.0 / 15.0;
    const double b43 = 32.0 / 9.0;
    const double b51 = 19372.0 / 6561.0;
    const double b52 = -25360.0 / 2187.0;
    const double b53 = 64448.0 / 6561.0;
    const double b54 = -212.0 / 729.0;
    const double b61 = 9017.0 / 3168.0;
    const double b62 = -355.0 / 33.0;
    const double b63 = 46732.0 / 5247.0;
    const double b64 = 49.0 / 176.0;
    const double b65 = -5103.0 / 18656.0;
    const double b71 = 35.0 / 384.0;
    const double b72 = 0.0;
    const double b73 = 500.0 / 1113.0;
    const double b74 = 125.0 / 192.0;
    const double b75 = -2187.0 / 6784.0;
    const double b76 = 11.0 / 84.0;

    const double c1 = 35.0 / 384.0;
    const double c2 = 0.0;
    const double c3 = 500.0 / 1113.0;
    const double c4 = 125.0 / 192.0;
    const double c5 = -2187.0 / 6784.0;
    const double c6 = 11.0 / 84.0;
    const double c7 = 0.0;

    const double d1 = 5179.0 / 57600.0;
    const double d2 = 0.0;
    const double d3 = 7571.0 / 16695.0;
    const double d4 = 393.0 / 640.0;
    const double d5 = -92097.0 / 339200.0;
    const double d6 = 187.0 / 2100.0;
    const double d7 = 1.0 / 40.0;

    // вычисление k1
    brusselator(t, z, k1);

    // вычисление k2
    for (size_t i = 0; i < z.size(); ++i) {
        z_temp[i] = z[i] + h * b21 * k1[i];
    }
    brusselator(t + a2 * h, z_temp, k2);

    // вычисление k3
    for (size_t i = 0; i < z.size(); ++i) {
        z_temp[i] = z[i] + h * (b31 * k1[i] + b32 * k2[i]);
    }
    brusselator(t + a3 * h, z_temp, k3);

    // вычисление k4
    for (size_t i = 0; i < z.size(); ++i) {
        z_temp[i] = z[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
    }
    brusselator(t + a4 * h, z_temp, k4);

    // вычисление k5
    for (size_t i = 0; i < z.size(); ++i) {
        z_temp[i] = z[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
    }
    brusselator(t + a5 * h, z_temp, k5);

    // вычисление k6
    for (size_t i = 0; i < z.size(); ++i) {
        z_temp[i] =
            z[i] + h * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
    }
    brusselator(t + a6 * h, z_temp, k6);

    // вычисление k7
    for (size_t i = 0; i < z.size(); ++i) {
        z_temp[i] = z[i] + h * (b71 * k1[i] + b72 * k2[i] + b73 * k3[i] + b74 * k4[i] +
                                 b75 * k5[i] + b76 * k6[i]);
    }
    brusselator(t + a7 * h, z_temp, k7);

    // вычисление ошибки
    double error = 0.0;
    for (size_t i = 0; i < z.size(); ++i) {
        double delta = h * ((d1 * k1[i]) + (d2 * k2[i]) + (d3 * k3[i]) + (d4 * k4[i]) +
                             (d5 * k5[i]) + (d6 * k6[i]) + (d7 * k7[i]));
        error += delta * delta;
    }
    error = sqrt(error / z.size());

    // выбор шага
    if (error > tolerance) {
        h *= 0.9 * pow(tolerance / error, 0.2);
    } else {
        h *= 1.1 * pow(tolerance / error, 0.2);
    }

    // обновление значений z
    for (size_t i = 0; i < z.size(); ++i) {
        z[i] += h * ((c1 * k1[i]) + (c2 * k2[i]) + (c3 * k3[i]) + (c4 * k4[i]) + (c5 * k5[i]) +
                      (c6 * k6[i]) + (c7 * k7[i]));
    }
}

void writeErrorToFile1(const vector<double> &z1, const vector<double> &z2, double time,
                       ofstream &file) {
    double error_x = fabs(z1[0] - z2[0]);
    double error_y = fabs(z1[1] - z2[1]);

    file << time << "\t" << error_x << "\t" << error_y << "\n";
}

void writeErrorToFile2(const vector<double> &z1, const vector<double> &z2, double time,
                       ofstream &file) {
    double error_x = fabs(z1[0] - z2[0]) / fabs(z2[0]);
    double error_y = fabs(z1[1] - z2[1]) / fabs(z2[1]);

    file << time << "\t" << error_x << "\t" << error_y << "\n";
}

void calculateTotalError(const string &filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file!\n";
        return;
    }

    double total_error_x = 0.0, total_error_y = 0.0;
    double t, error_x, error_y;

    while (file >> t >> error_x >> error_y) {
        total_error_x += error_x;
        total_error_y += error_y;
    }

    file.close();

    std::cout << "Total error for file " << filename << ":\n";
    std::cout << "Error on x: " << total_error_x << "\n";
    std::cout << "Error on y: " << total_error_y << "\n";
}

int main() {
    // начальные условия
    vector<double> z1 = {1.0, 1.0}; // x(0) = 1.0, y(0) = 1.0
    vector<double> z2 = {1.0, 1.0}; // x(0) = 1.0, y(0) = 1.0

    // временные параметры
    double t = 0.0;
    double t_end = 20.0;
    double h = 0.01;
    double tolerance = 1e-4;

    // открываем файл для записи результатов
    ofstream outputFile_rk4("brusselator_output.txt");
    ofstream errorFile1("error1.txt");
    ofstream errorFile2("error2.txt");
    if (!outputFile_rk4.is_open() || !errorFile1.is_open() || !errorFile2.is_open()) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    // заголовок файла
    outputFile_rk4 << "Time\tx\ty\n";

    // численное решение системы
    while (t <= t_end) {
        outputFile_rk4 << t << "\t" << z1[0] << "\t" << z1[1] << "\n";
        rungeKutta4(t, h, z1);
        dormandPrince5(t, h, z2, tolerance);

        writeErrorToFile1(z1, z2, t, errorFile1);
        writeErrorToFile2(z1, z2, t, errorFile2);

        t += h;
    }

    outputFile_rk4.close();
    errorFile1.close();
    errorFile2.close();
    cout << "The results are written to a file brusselator_output.txt\n\n";
    cout << "File error1.txt: absolute errors\n";
    calculateTotalError("error1.txt");
    cout << "\n\nFile error2.txt: relative errors\n";
    calculateTotalError("error2.txt");
    return 0;
}