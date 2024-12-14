#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;

int nextPowerOfTwo(int N) {
    int power = 1;
    while (power < N) {
        power <<= 1;
    }
    return power;
}

void createMatrix(int K, double**& matrix) {
    matrix = new double*[K];
    for (int i = 0; i < K; ++i) {
        matrix[i] = new double[K]();
    }
}

void fillMatrixRandomly(int N, double** matrix) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            matrix[i][j] = rand() % 101;
        }
    }
}

void addMatrices(int K, double** A, double** B, double** C) {
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

void dedMatrices(int K, double** A, double** B, double** C) {
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

void multiplyMatricesStandard(int N, double** A, double** B, double** C) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < N; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void strassen(int K, double** A, double** B, double** C) {
    if (K == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }

    int newK = K / 2;
    double** A11 = new double*[newK];
    double** A12 = new double*[newK];
    double** A21 = new double*[newK];
    double** A22 = new double*[newK];
    double** B11 = new double*[newK];
    double** B12 = new double*[newK];
    double** B21 = new double*[newK];
    double** B22 = new double*[newK];
    double** C11 = new double*[newK];
    double** C12 = new double*[newK];
    double** C21 = new double*[newK];
    double** C22 = new double*[newK];
    double** P1 = new double*[newK];
    double** P2 = new double*[newK];
    double** P3 = new double*[newK];
    double** P4 = new double*[newK];
    double** P5 = new double*[newK];
    double** P6 = new double*[newK];
    double** P7 = new double*[newK];
    double** D1 = new double*[newK];
    double** D2 = new double*[newK];
    for (int i = 0; i < newK; ++i) {
        A11[i] = new double[newK]();
        A12[i] = new double[newK]();
        A21[i] = new double[newK]();
        A22[i] = new double[newK]();
        B11[i] = new double[newK]();
        B12[i] = new double[newK]();
        B21[i] = new double[newK]();
        B22[i] = new double[newK]();
        C11[i] = new double[newK]();
        C12[i] = new double[newK]();
        C21[i] = new double[newK]();
        C22[i] = new double[newK]();
        P1[i] = new double[newK]();
        P2[i] = new double[newK]();
        P3[i] = new double[newK]();
        P4[i] = new double[newK]();
        P5[i] = new double[newK]();
        P6[i] = new double[newK]();
        P7[i] = new double[newK]();
        D1[i] = new double[newK]();
        D2[i] = new double[newK]();
    }

    for (int i = 0; i < newK; ++i) {
        for (int j = 0; j < newK; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + newK];
            A21[i][j] = A[i + newK][j];
            A22[i][j] = A[i + newK][j + newK];
            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + newK];
            B21[i][j] = B[i + newK][j];
            B22[i][j] = B[i + newK][j + newK];
        }
    }

    addMatrices(newK, A11, A22, D1);
    addMatrices(newK, B11, B22, D2);
    strassen(newK, D1, D2, P1);

    dedMatrices(newK, A12, A22, D1);
    addMatrices(newK, B21, B22, D2);
    strassen(newK, D1, D2, P2);

    dedMatrices(newK, A21, A11, D1);
    addMatrices(newK, B11, B12, D2);
    strassen(newK, D1, D2, P3);

    addMatrices(newK, A11, A12, D1);
    strassen(newK, D1, B22, P4);

    addMatrices(newK, A21, A22, D1);
    strassen(newK, D1, B11, P5);

    dedMatrices(newK, B21, B11, D1);
    strassen(newK, A22, D1, P6);

    dedMatrices(newK, B12, B22, D1);
    strassen(newK, A11, D1, P7);

    addMatrices(newK, P1, P2, D1);
    dedMatrices(newK, P6, P4, D2);
    addMatrices(newK, D1, D2, C11);

    addMatrices(newK, P1, P3, D1);
    dedMatrices(newK, P7, P5, D2);
    addMatrices(newK, D1, D2, C22);

    addMatrices(newK, P7, P4, C12);
    addMatrices(newK, P6, P5, C21);

    for (int i = 0; i < newK; ++i) {
        for (int j = 0; j < newK; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + newK] = C12[i][j];
            C[i + newK][j] = C21[i][j];
            C[i + newK][j + newK] = C22[i][j];
        }
    }

    for (int i = 0; i < newK; ++i) {
        delete[] A11[i];
        delete[] A12[i];
        delete[] A21[i];
        delete[] A22[i];
        delete[] B11[i];
        delete[] B12[i];
        delete[] B21[i];
        delete[] B22[i];
        delete[] C11[i];
        delete[] C12[i];
        delete[] C21[i];
        delete[] C22[i];
        delete[] P1[i];
        delete[] P2[i];
        delete[] P3[i];
        delete[] P4[i];
        delete[] P5[i];
        delete[] P6[i];
        delete[] P7[i];
        delete[] D1[i];
        delete[] D2[i];
    }
    delete[] A11;
    delete[] A12;
    delete[] A21;
    delete[] A22;
    delete[] B11;
    delete[] B12;
    delete[] B21;
    delete[] B22;
    delete[] C11;
    delete[] C12;
    delete[] C21;
    delete[] C22;
    delete[] P1;
    delete[] P2;
    delete[] P3;
    delete[] P4;
    delete[] P5;
    delete[] P6;
    delete[] P7;
    delete[] D1;
    delete[] D2;
}

int main() {
    srand(static_cast<unsigned>(time(0)));

    int N;
    cout << "Enter the size of the matrix N: ";
    cin >> N;

    int K = nextPowerOfTwo(N);
    double** A = nullptr;
    double** B = nullptr;
    double** C = nullptr;
    double** C_standard = nullptr;

    createMatrix(K, A);
    createMatrix(K, B);
    createMatrix(K, C);
    createMatrix(K, C_standard);

    fillMatrixRandomly(N, A);
    fillMatrixRandomly(N, B);

    cout << "Matrix A:" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << setw(5) << A[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Matrix B:" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << setw(5) << B[i][j] << " ";
        }
        cout << endl;
    }

    strassen(K, A, B, C);
    multiplyMatricesStandard(N, A, B, C_standard);

    cout << "Resulting matrix C (Strassen):" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << setw(10) << C[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Resulting matrix C (Standard Multiplication):" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << setw(10) << C_standard[i][j] << " ";
        }
        cout << endl;
    }

    for (int i = 0; i < K; ++i) {
        delete[] A[i];
        delete[] B[i];
        delete[] C[i];
        delete[] C_standard[i];
    }
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] C_standard;

    return 0;
}
