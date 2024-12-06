#include <vector>
#include <cmath>
#include <iostream>


// this entire code section was written with chat gpt

//  I used the following prompt with gpt 4-o
//  write c++ code to solve the upcomming matrix calculations using conjugate gradient method: 
//  Q * x + dx = 0; where Q is a n x n matrix, x and dx are of length n. 
//  the values of Q & dx are already known, find x 


// Function to perform matrix-vector multiplication O(vec.size()^2)
std::vector<double> multiplyMatrixVector(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector) {
    int size = vector.size();
    std::vector<double> result(size, 0.0);
    for (int row = 0; row < size; ++row) {
        for (int col = 0; col < size; ++col) {
            result[row] += matrix[row][col] * vector[col];
        }
    }
    return result;
}

// Function to compute the dot product of two vectors (O(vec1.size()))
double computeDotProduct(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    double result = 0.0;
    for (size_t index = 0; index < vector1.size(); ++index) {
        result += vector1[index] * vector2[index];
    }
    return result;
}

// Function to add two vectors (O(vec1.size()))
std::vector<double> addVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    int size = vector1.size();
    std::vector<double> result(size);
    for (int index = 0; index < size; ++index) {
        result[index] = vector1[index] + vector2[index];
    }
    return result;
}

// Function to subtract two vectors (O(vec1.size()))
std::vector<double> subtractVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    int size = vector1.size();
    std::vector<double> result(size);
    for (int index = 0; index < size; ++index) {
        result[index] = vector1[index] - vector2[index];
    }
    return result;
}

// Function to scale a vector by a scalar O(vec.size())
std::vector<double> scaleVector(const std::vector<double>& vector, double scalar) {
    int size = vector.size();
    std::vector<double> result(size);
    for (int index = 0; index < size; ++index) {
        result[index] = vector[index] * scalar;
    }
    return result;
}

// Conjugate Gradient method to solve Q * x + dx = 0
// Conjugate Gradient method to solve A * solution + rhs = 0
std::vector<double> solveConjugateGradient(const std::vector<std::vector<double>>& matrixA, const std::vector<double>& vectorRHS, double tolerance = 1e-6, int maxIterations = 1000) {
    int size = vectorRHS.size();
    std::vector<double> solution(size, 0.0); // Initialize solution vector
    std::vector<double> residual = scaleVector(vectorRHS, -1.0); // Compute initial residual
    std::vector<double> direction = residual;
    double residualSquaredOld = computeDotProduct(residual, residual);

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        std::cout << "Running Iteration " << iteration << " out of " << maxIterations << " Total Iterations" << std::endl;
        std::vector<double> matrixDirection = multiplyMatrixVector(matrixA, direction);
        double stepSize = residualSquaredOld / computeDotProduct(direction, matrixDirection);
        solution = addVectors(solution, scaleVector(direction, stepSize));
        residual = subtractVectors(residual, scaleVector(matrixDirection, stepSize));
        double residualSquaredNew = computeDotProduct(residual, residual);

        if (sqrt(residualSquaredNew) < tolerance) {
            std::cout << "Finished on Iteration: " << iteration << std::endl;
            break;
        }

        direction = addVectors(residual, scaleVector(direction, residualSquaredNew / residualSquaredOld));
        residualSquaredOld = residualSquaredNew;
    }

    return solution;
}