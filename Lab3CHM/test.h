#pragma once
#include "solve.h"

using matrix = std::vector<std::vector<double>>;

struct TestData
{
    double absMinLambda;
    int minLambdaNum;
    std::vector<double> vecX;
    OutParams out;
};

const double e = 1.0e-17;
const double q = 1.0e10;

matrix multiplyMatrices(const matrix& matrix1, const matrix& matrix2)
{
    if (matrix1[0].size() != matrix2.size()) {
        return matrix();
    }

    matrix result(matrix1.size(), std::vector<double>(matrix2[0].size(), 0));

    for (size_t i = 0; i < matrix1.size(); ++i) {
        for (size_t j = 0; j < matrix2[0].size(); ++j) {
            for (size_t k = 0; k < matrix2.size(); ++k) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

matrix multiplyMatrixByScalar(const matrix& inMatrix, double scalar)
{
    size_t rows = inMatrix.size();
    size_t cols = inMatrix[0].size();

    matrix result(rows, std::vector<double>(cols, 0));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[i][j] = inMatrix[i][j] * scalar;
        }
    }

    return result;
}

matrix subtractMatrices(const matrix& matrix1, const matrix& matrix2) {
    if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
        return matrix();
    }

    matrix result(matrix1.size(), std::vector<double>(matrix1[0].size(), 0));

    for (size_t i = 0; i < matrix1.size(); ++i) {
        for (size_t j = 0; j < matrix1[0].size(); ++j) {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }

    return result;
}

matrix getLambdai(int min, int max, int size, TestData& tInput)
{
    tInput.absMinLambda = max;

    matrix result(size, std::vector<double>(size, 0));
    for (int i = 0; i < size; i++)
    {
        result[i][i] = (double)(min * 1000 + rand() % (int)((max - min) * 1000)) / 1000.0;
        if (abs(result[i][i]) < tInput.absMinLambda)
        {
            tInput.absMinLambda = abs(result[i][i]);
            tInput.minLambdaNum = i;
        }
    }
    return result;
}

matrix multiplyVectorsToMat(const std::vector<double>& vector1, const std::vector<double>& vector2)
{
    matrix result(vector1.size(), std::vector<double>(vector2.size(), 0));

    for (size_t i = 0; i < vector1.size(); ++i) {
        for (size_t j = 0; j < vector2.size(); ++j) {
            result[i][j] = vector1[i] * vector2[j];
        }
    }

    return result;
}

std::vector<double> getNormedW(int min, int max, int size)
{
    std::vector<double> result(size, 0.0);
    for (int i = 0; i < size; i++)
        result[i] = (double)(min * 1000 + rand() % (int)((max - min) * 1000)) / 1000.0;

    double norm = 0.0;
    for (int i = 0; i < size; i++)
        norm += (result[i] * result[i]);
    norm = sqrt(norm);

    for (int i = 0; i < size; i++)
        result[i] /= norm;
    return result;
}

matrix transposeMatrix(const matrix& inMatrix)
{
    size_t rows = inMatrix.size();
    size_t cols = inMatrix[0].size();

    matrix result(cols, std::vector<double>(rows, 0));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[j][i] = inMatrix[i][j];
        }
    }

    return result;
}

matrix getE(int size)
{
    matrix result(size, std::vector<double>(size, 0.0));
    for (int i = 0; i < size; i++)
    {
        result[i][i] = 1.0;
    }
    return result;
}

matrix getHmatrix(int size, int min, int max)
{
    std::vector<double> w = getNormedW(min, max, size);
    matrix wwt2 = multiplyMatrixByScalar(multiplyVectorsToMat(w, w), 2);
    return subtractMatrices(getE(size), wwt2);
}

double calcLambdaError(double lambda, double calcLambda)
{
    double tmp = abs(calcLambda) - abs(lambda);
    return (lambda < e) ? abs(tmp) : abs(tmp / lambda);
}

double calcXError(const std::vector<double>& x, const std::vector<double>& calcX)
{
    double error = 0.0;
    for (int i = 0; i < x.size(); i++)
    {
        error = std::max(error, calcLambdaError(x[i], calcX[i]));
    }
    return error;
}

std::vector<double> colToVec(matrix mat, int colNum)
{
    std::vector<double> result(mat.size(), 0.0);
    for (int i = 0; i < mat.size(); i++)
    {
        result[i] = mat[i][colNum];
    }
    return result;
}

TestData testIteration(int min, int max, double eps, int size)
{
    matrix hmat = getHmatrix(size, min, max);
    //matrix hmat = GenerateHouseholderMatrix(getNormedW(min, max, size));
    TestData data;
    matrix lambdas = getLambdai(min, max, size, data);
    data.vecX = colToVec(hmat, data.minLambdaNum);

    matrix A = multiplyMatrices(multiplyMatrices(hmat, lambdas), transposeMatrix(hmat));

    InParams in;
    in.N = size;
    in.lambda_epsilon = in.g_epsilon = eps;
    in.A = A;
    in.M = 10000;

    OutParams out;
    Solve(in, out);
    data.out = out;

    return data;
}

void test(int testIters, int min, int max, double eps, int size)
{
    int avgItersCount = 0;
    double avgr = 0.0;
    double avgLambdaError = 0.0;
    double avgXError = 0.0;
    for (int i = 0; i < testIters; i++)
    {
        TestData data = testIteration(min, max, eps, size);
        avgItersCount += data.out.K;
        avgr += data.out.r;
        avgLambdaError += calcLambdaError(data.absMinLambda, data.out.lambda);
        avgXError += calcXError(data.vecX, data.out.x);
    }
    avgItersCount /= testIters;
    avgr /= testIters;
    avgLambdaError /= testIters;
    avgXError /= testIters;

    std::cout << "N: " << size << "  lambda:" << min << ":" << max << "  eps:" << eps
        << "   Ср собственных:" << avgLambdaError << "  ср векторов:" << avgXError << "   avgR:" << avgr << "   avgIters:" << avgItersCount << '\n';
}

void testsRun()
{
    test(15, -2, 2, 1e-5, 10);
    test(15, -2, 2, 1e-8, 10);
    test(15, -50, 50, 1e-5, 10);
    test(15, -50, 50, 1e-8, 10);

    test(15, -2, 2, 1e-5, 30);
    test(15, -2, 2, 1e-8, 30);
    test(15, -50, 50, 1e-5, 30);
    test(15, -50, 50, 1e-8, 30);

    test(15, -2, 2, 1e-5, 50);
    test(15, -2, 2, 1e-8, 50);
    test(15, -50, 50, 1e-5, 50);
    test(15, -50, 50, 1e-8, 50);
}