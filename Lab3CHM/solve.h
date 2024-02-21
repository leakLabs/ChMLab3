#pragma once

#include <iostream>
#include <fstream>
#include <vector>

struct InParams {
    int N;
    std::vector<std::vector<double>> A;
    double lambda_epsilon;
    double g_epsilon;
    int M;
};

struct OutParams {
    int IER;
    double lambda;
    std::vector<double> x;
    int K;
    double r;
};

// Метод Гаусса для решения системы линейных уравнений Ax = f
std::vector<double> gaussElimination(const std::vector<std::vector<double>>& A, const std::vector<double>& f) {
    int n = A.size();

    // Создаем расширенную матрицу [A | f]
    std::vector<std::vector<double>> augmentedMatrix(n, std::vector<double>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmentedMatrix[i][j] = A[i][j];
        }
        augmentedMatrix[i][n] = f[i];
    }

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; ++i) {
        double pivot = augmentedMatrix[i][i];
        for (int j = i; j < n + 1; ++j) {
            augmentedMatrix[i][j] /= pivot;
        }

        // Вычитание i-й строки из оставшихся строк
        for (int k = i + 1; k < n; ++k) {
            double factor = augmentedMatrix[k][i];
            for (int j = i; j < n + 1; ++j) {
                augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
            }
        }
    }

    // Обратный ход метода Гаусса
    std::vector<double> solution(n);
    for (int i = n - 1; i >= 0; --i) {
        solution[i] = augmentedMatrix[i][n];
        for (int j = i + 1; j < n; ++j) {
            solution[i] -= augmentedMatrix[i][j] * solution[j];
        }
    }

    return solution;
}

double calculate_error(const std::vector<std::vector<double>>& A, const std::vector<double>& x, double lambda) {
    std::vector<double> result;

    // Перемножение матрицы A и вектора x и вычитание lambdax
    for (size_t i = 0; i < A.size(); ++i) {
        double element = 0.0;
        for (size_t j = 0; j < A[i].size(); ++j) {
            element += A[i][j] * x[j];
        }
        result.push_back(element - lambda * x[i]);
    }

    // Вычисление погрешности по первой норме
    double max_abs_error = 0.0;
    for (size_t i = 0; i < result.size(); ++i) {
        max_abs_error = std::max(max_abs_error, std::abs(result[i]));
    }

    return max_abs_error;
}

//Нормализация вектора
std::vector<double> Normalize(const std::vector<double>& vec) {
    if (vec.size() == 0) {
        std::cerr << "Error: Cannot normalize a zero-length vector." << std::endl;
        return {};
    }

    std::vector<double> normalizedVec = vec;
    double sumOfSquares = 0.0;

    for (double val : normalizedVec) {
        sumOfSquares += val * val;
    }

    double norm = std::sqrt(sumOfSquares);

    for (double& val : normalizedVec) {
        val /= norm;
    }

    return normalizedVec;
}

//Генерация матрицы хаусхолдера по вектору w
std::vector<std::vector<double>> GenerateHouseholderMatrix(const std::vector<double>& w) {
    int n = w.size();

    if (n == 0) {
        std::cerr << "Error: a zero-length vector" << std::endl;
        return {};
    }

    // Нормализация вектора w
    std::vector<double> normalizedW = Normalize(w);

    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        H[i][i] = 1.0;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            H[i][j] -= 2.0 * normalizedW[i] * normalizedW[j];
        }
    }

    return H;
}

//Создать диагональную матрицу из лямбда вектора
std::vector<std::vector<double>> CreateDiagonalMatrix(const std::vector<double>& lambdas) {
    int n = lambdas.size();

    if (n == 0) {
        std::cerr << "Error: a zero-length vector" << std::endl;
        return {};
    }

    std::vector<std::vector<double>> diagonalMatrix(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        diagonalMatrix[i][i] = lambdas[i];
    }

    return diagonalMatrix;
}


// Вычисление A = HMH^T
std::vector<std::vector<double>> TransformMatrix(const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& H) {
    int n = M.size();


    std::vector<std::vector<double>> res1(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                res1[i][j] += H[i][k] * M[k][j];
            }
        }
    }
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {

            for (int k = 0; k < n; ++k) {
                result[i][j] += res1[i][k] * H[j][k];
            }
        }
    }

    return result;

}

std::vector<std::vector<double>> CreateA(const std::vector<double>& w, const std::vector<double> lambdas, std::vector<std::vector<double>>& H) {
    std::vector<std::vector<double>> A;
    H = GenerateHouseholderMatrix(w);
    A = TransformMatrix(CreateDiagonalMatrix(lambdas), H);
    return A;
}


void ReadInputFromFile(const std::string& filename, InParams& params) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error while open file.\n";
        return;
    }

    inputFile >> params.N;

    params.A.resize(params.N, std::vector<double>(params.N, 0.0));
    for (int i = 0; i < params.N; ++i) {
        for (int j = 0; j < params.N; ++j) {
            inputFile >> params.A[i][j];
        }
    }

    inputFile >> params.lambda_epsilon;
    inputFile >> params.g_epsilon;
    inputFile >> params.M;

    inputFile.close();
}

void print(OutParams& params) {
    std::cout << "Собственное значение: " << params.lambda << "\nсобственный вектор: ";
    for (double num : params.x)
        std::cout << num << ' ';
    std::cout << "\nЧисло итераций: " << params.K << "   Мера точности :" << params.r << '\n';
}

//Вычисление угла между углами
double Angle(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    double dotProduct = 0.0;
    double magnitudeVector1 = 0.0;
    double magnitudeVector2 = 0.0;

    for (size_t i = 0; i < vector1.size(); ++i) {
        dotProduct += -vector1[i] * vector2[i];
        magnitudeVector1 += vector1[i] * vector1[i];
        magnitudeVector2 += vector2[i] * vector2[i];
    }

    return std::acos(dotProduct / (sqrt(magnitudeVector1) * sqrt(magnitudeVector2)));
}

double GetSigma(const std::vector<double>& v, const std::vector<double>& x) {
    int size = v.size();
    double result = 0;
    for (int i = 0; i < size; ++i) {
        result += v[i] * x[i];
    }
    return  1 / result;
}

void Solve(const InParams& in, OutParams& out) {
    out.K = 3;

    std::vector<double> prevX(in.A.size(), 1);
    std::vector<double> prevV = Normalize(prevX);
    double prevSigma;

    std::vector<double> currX = gaussElimination(in.A, prevV);
    std::vector<double> currV = Normalize(currX);
    double currSigma = GetSigma(prevV, currX);

    prevX = currX;
    prevV = currV;
    prevSigma = currSigma;
    currX = gaussElimination(in.A, prevV);
    currV = Normalize(currX);
    currSigma = GetSigma(prevV, currX);

    while ((Angle(currX, prevX) > in.g_epsilon || std::abs(currSigma - prevSigma) > in.lambda_epsilon) && out.K < in.M) {
        double a = Angle(currX, prevX);
        prevX = currX;
        prevV = currV;
        prevSigma = currSigma;
        currX = gaussElimination(in.A, prevV);
        currV = Normalize(currX);
        currSigma = GetSigma(prevV, currX);
        ++out.K;
    }

    if (out.K >= in.M) {
        out.IER = -1;
    }
    else {
        out.IER = 0;
    }
    out.lambda = currSigma;
    out.x = currV;

    if (out.K % 2 == 0) {
        for (double& val : out.x) {
            val *= -1;
        }
    }

    out.r = calculate_error(in.A, out.x, out.lambda);
}