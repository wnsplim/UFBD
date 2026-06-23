#include "GtrRateModel.hpp"
#include "Eigen/Dense"

#include <cmath>

GtrRateModel::GtrRateModel(int n) :
    numStates(n),
    eigenvalue(n),
    cijk(n * n * n)
{
}

void GtrRateModel::setParameters(const std::vector<double>& exch, const std::vector<double>& freq){
    int n = numStates;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(n, n);
    int k = 0;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            Q(i, j) = exch[k] * freq[j];
            Q(j, i) = exch[k] * freq[i];
            k++;
        }
    }
    for(int i = 0; i < n; i++){
        double s = 0.0;
        for(int j = 0; j < n; j++)
            if(j != i) s += Q(i, j);
        Q(i, i) = -s;
    }
    double mu = 0.0;
    for(int i = 0; i < n; i++)
        mu += freq[i] * -Q(i, i);
    Q /= mu;

    Eigen::MatrixXd B(n, n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            B(i, j) = std::sqrt(freq[i]) * Q(i, j) / std::sqrt(freq[j]);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
    const Eigen::VectorXd& lam = es.eigenvalues();
    const Eigen::MatrixXd& U = es.eigenvectors();
    for(int s = 0; s < n; s++)
        eigenvalue[s] = lam(s);

    double* p = &cijk[0];
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            double c = std::sqrt(freq[j] / freq[i]);
            for(int s = 0; s < n; s++)
                *(p++) = U(i, s) * U(j, s) * c;
        }
    }
}

void GtrRateModel::transitionProbabilities(double t, double* P) const {
    int n = numStates;
    std::vector<double> ev(n);
    for(int s = 0; s < n; s++)
        ev[s] = std::exp(eigenvalue[s] * t);
    const double* c = &cijk[0];
    for(int i = 0; i < n; i++){
        double rowsum = 0.0;
        for(int j = 0; j < n; j++){
            double sum = 0.0;
            for(int s = 0; s < n; s++)
                sum += *(c++) * ev[s];
            if(sum < 0.0) sum = 0.0;
            P[i * n + j] = sum;
            rowsum += sum;
        }
        for(int j = 0; j < n; j++)
            P[i * n + j] /= rowsum;
    }
}
