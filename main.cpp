#include <iostream>
#include "vector"
#include "cmath"
#include "fstream"
#include "sstream"
#include "string"

void w_out(std::vector<std::vector<double>> const& M, std::ofstream& out){
    for (int i = 0; i < M.size(); ++i) {
        out << "{ ";
        for (int j = 0; j < M.size(); j++) {
            if (M[i][j] == 0) {
                out << 0;
            } else if (M[i][j] == 1){
                out << 1;
            } else {
                std::stringstream s;
                s << M[i][j];
                std::string str = s.str();
                if (str.find('e',0) != -1) {
                    str.replace(str.find('e', 0), 1, "*10^");
                    if (str.find('+',0) != -1){
                        str.replace(str.find('+', 0), 1, "");
                    }
                    out << str;
                } else {
                    out << M[i][j];
                }
            }
            if (j != M.size()-1){
                out << ", ";
            }
        }
        out << " } ";
        if (i != M.size()-1){
            out << std::endl;
        }
    }
}

void v_out(std::vector<double> const& V, std::ofstream& out){
    for (int i = 0; i < V.size(); ++i) {
        if (V[i] == 0) {
            out << 0 << "  ";
        } else if (V[i] == 1){
            out << 1 << "  ";
        } else {
            std::stringstream s;
            s << V[i];
            std::string str = s.str();
            if (str.find('e',0) != -1) {
                str.replace(str.find('e', 0), 1, "*10^");
                if (str.find('+',0) != -1){
                    str.replace(str.find('+', 0), 1, "");
                }
                out << str << "  ";
            } else {
                out << V[i] << "  ";
            }
        }
        if (i != V.size()-1){
            out << std::endl;
        }
    }
}

int main() {
    double nu = 0.3;
    double a = 0.1;
    double b = 0.2;
    double U = 0.00001;
    double E = 1.1*std::pow(10,10);
    int n = 71;
    int m = 71;
    double h = (b-a)/(n-1);
    double tau = 2*M_PI/(m-1);
    std::vector<double> B(n*m*2);
    for (double & i : B) {
        i = 0;
    }
    std::vector<std::vector<double>> M(n*m*2);
    for (auto & i : M) {
        i = std::vector<double>(n*m*2);
        for (int j = 0; j < n*m*2; ++j) {
            i[j] = 0;
        }
    }
    M[0][0] = -E*a*tau/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[0][m] = E*h*nu*tau/(h*a*tau-h*a*nu*nu*tau);
    M[0][2*m] = E*a*tau/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[0][n*m+2*m-1] = -E*h*nu/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[0][n*m+m+1] = E*h*nu/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[1][0] = -E*a*tau*nu/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[1][m] = E*h*tau/(h*a*tau-h*a*nu*nu*tau);
    M[1][2*m] = E*a*nu*tau/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[1][n*m+2*m-1] = -E*h/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[1][n*m+m+1] = E*h/(2*h*a*tau-2*h*a*nu*nu*tau);
    for(int j = 1; j < m-1; ++j){
        M[2*j][j] = -E*a*tau/(2*h*a*tau-2*h*a*nu*nu*tau);
        M[2*j][m+j] = E*h*nu*tau/(h*a*tau-h*a*nu*nu*tau);
        M[2*j][2*m+j] = E*a*tau/(2*h*a*tau-2*h*a*nu*nu*tau);
        M[2*j][n*m+m+j-1] = -E*h*nu/(2*h*a*tau-2*h*a*nu*nu*tau);
        M[2*j][n*m+m+j+1] = E*h*nu/(2*h*a*tau-2*h*a*nu*nu*tau);
        M[2*j+1][j] = -E*a*tau*nu/(2*h*a*tau-2*h*a*nu*nu*tau);
        M[2*j+1][m+j] = E*h*tau/(h*a*tau-h*a*nu*nu*tau);
        M[2*j+1][2*m+j] = E*a*nu*tau/(2*h*a*tau-2*h*a*nu*nu*tau);
        M[2*j+1][n*m+m+j-1] = -E*h/(2*h*a*tau-2*h*a*nu*nu*tau);
        M[2*j+1][n*m+m+j+1] = E*h/(2*h*a*tau-2*h*a*nu*nu*tau);
    }
    M[2*(m-1)][m-1] = -E*a*tau/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[2*(m-1)][2*m-1] = E*h*nu*tau/(h*a*tau-h*a*nu*nu*tau);
    M[2*(m-1)][3*m-1] = E*a*tau/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[2*(m-1)][n*m+m] = -E*h*nu/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[2*(m-1)][n*m+2*m+1] = E*h*nu/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[2*(m-1)+1][m-1] = -E*a*tau*nu/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[2*(m-1)+1][2*m-1] = E*h*tau/(h*a*tau-h*a*nu*nu*tau);
    M[2*(m-1)+1][3*m-1] = E*a*nu*tau/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[2*(m-1)+1][n*m+m] = -E*h/(2*h*a*tau-2*h*a*nu*nu*tau);
    M[2*(m-1)+1][n*m+2*m+1] = E*h/(2*h*a*tau-2*h*a*nu*nu*tau);
    for (int i = 1; i < n-1; ++i) {
        for (int j = 1; j < m-1; ++j) {
            M[2*m*i+2*j][(i-1)*m+j-1] = 0;
            M[2*m*i+2*j][(i-1)*m+j] = 2*a*a/(h*h) - (a/h) - i +4*a*i/h +2*i*i;
            M[2*m*i+2*j][(i-1)*m+j+1] = 0;
            M[2*m*i+2*j][i*m+j-1] = 1/(tau*tau) - (nu/(tau*tau));
            M[2*m*i+2*j][i*m+j] = -2-4*a*a/(h*h) - (8*a*i)/h - 4*i*i - 2/(tau*tau)+2*nu/(tau*tau);
            M[2*m*i+2*j][i*m+j+1] = 1/(tau*tau) - (nu/(tau*tau));
            M[2*m*i+2*j][(i+1)*m+j-1] = 0;
            M[2*m*i+2*j][(i+1)*m+j] = 2*a*a/(h*h) - (a/h) - i +4*a*i/h +2*i*i;
            M[2*m*i+2*j][(i+1)*m+j+1] = 0;
            M[2*m*i+2*j][n*m+(i-1)*m+j-1] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
            M[2*m*i+2*j][n*m+(i-1)*m+j] = 0;
            M[2*m*i+2*j][n*m+(i-1)*m+j+1] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
            M[2*m*i+2*j][n*m+i*m+j-1] = 3/(2*tau)-nu/(2*tau);
            M[2*m*i+2*j][n*m+i*m+j] = 0;
            M[2*m*i+2*j][n*m+i*m+j+1] = -3/(2*tau)+nu/(2*tau);
            M[2*m*i+2*j][n*m+(i+1)*m+j-1] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
            M[2*m*i+2*j][n*m+(i+1)*m+j] = 0;
            M[2*m*i+2*j][n*m+(i+1)*m+j+1] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
            //
            M[2*m*i+2*j+1][(i-1)*m+j-1] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
            M[2*m*i+2*j+1][(i-1)*m+j] = 0;
            M[2*m*i+2*j+1][(i-1)*m+j+1] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
            M[2*m*i+2*j+1][i*m+j-1] = -3/(2*tau)+nu/(2*tau);
            M[2*m*i+2*j+1][i*m+j] = 0;
            M[2*m*i+2*j+1][i*m+j+1] = 3/(2*tau)-nu/(2*tau);
            M[2*m*i+2*j+1][(i+1)*m+j-1] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
            M[2*m*i+2*j+1][(i+1)*m+j] = 0;
            M[2*m*i+2*j+1][(i+1)*m+j+1] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
            M[2*m*i+2*j+1][n*m+(i-1)*m+j-1] = 0;
            M[2*m*i+2*j+1][n*m+(i-1)*m+j] = -(a+h*i)*(2*a*(nu-1)*nu+h*(-1+(-1+2*i*(nu-1))))/(2*h*h);
            M[2*m*i+2*j+1][n*m+(i-1)*m+j+1] = 0;
            M[2*m*i+2*j+1][n*m+i*m+j-1] = 2/(tau*tau);
            M[2*m*i+2*j+1][n*m+i*m+j] = -1+nu+2*(a+h*i)*(a+h*i)*(nu-1)*nu/(h*h)-4/(tau*tau);
            M[2*m*i+2*j+1][n*m+i*m+j+1] = 2/(tau*tau);;
            M[2*m*i+2*j+1][n*m+(i+1)*m+j-1] = 0;
            M[2*m*i+2*j+1][n*m+(i+1)*m+j] = -(a+h*i)*(2*a*(nu-1)*nu+h*(-1+(-1+2*i*(nu-1))*nu))/(2*h*h);
            M[2*m*i+2*j+1][n*m+(i+1)*m+j+1] = 0;
        }
    }
    for (int i = 1; i < n; ++i) {
        M[2*m*i][(i-1)*m+m-1] = 0;
        M[2*m*i][(i-1)*m] = 2*a*a/(h*h) - (a/h) - i +4*a*i/h +2*i*i;
        M[2*m*i][(i-1)*m+1] = 0;
        M[2*m*i][i*m+m-1] = 1/(tau*tau) - (nu/(tau*tau));
        M[2*m*i][i*m] = -2-4*a*a/(h*h) - (8*a*i)/h - 4*i*i - 2/(tau*tau)+2*nu/(tau*tau);
        M[2*m*i][i*m+1] = 1/(tau*tau) - (nu/(tau*tau));
        M[2*m*i][(i+1)*m+m-1] = 0;
        M[2*m*i][(i+1)*m] = 2*a*a/(h*h) - (a/h) - i +4*a*i/h +2*i*i;
        M[2*m*i][(i+1)*m+1] = 0;
        M[2*m*i][n*m+(i-1)*m+m-1] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
        M[2*m*i][n*m+(i-1)*m] = 0;
        M[2*m*i][n*m+(i-1)*m+1] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
        M[2*m*i][n*m+i*m+m-1] = 3/(2*tau)-nu/(2*tau);
        M[2*m*i][n*m+i*m] = 0;
        M[2*m*i][n*m+i*m+1] = -3/(2*tau)+nu/(2*tau);
        M[2*m*i][n*m+(i+1)*m+m-1] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
        M[2*m*i][n*m+(i+1)*m] = 0;
        M[2*m*i][n*m+(i+1)*m+1] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
        //
        M[2*m*i+1][(i-1)*m+m-1] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
        M[2*m*i+1][(i-1)*m] = 0;
        M[2*m*i+1][(i-1)*m+1] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
        M[2*m*i+1][i*m+m-1] = -3/(2*tau)+nu/(2*tau);
        M[2*m*i+1][i*m] = 0;
        M[2*m*i+1][i*m+1] = 3/(2*tau)-nu/(2*tau);
        M[2*m*i+1][(i+1)*m+m-1] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
        M[2*m*i+1][(i+1)*m] = 0;
        M[2*m*i+1][(i+1)*m+1] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
        M[2*m*i+1][n*m+(i-1)*m+m-1] = 0;
        M[2*m*i+1][n*m+(i-1)*m] = -(a+h*i)*(2*a*(nu-1)*nu+h*(-1+(-1+2*i*(nu-1))))/(2*h*h);
        M[2*m*i+1][n*m+(i-1)*m+1] = 0;
        M[2*m*i+1][n*m+i*m+m-1] = 2/(tau*tau);
        M[2*m*i+1][n*m+i*m] = -1+nu+2*(a+h*i)*(a+h*i)*(nu-1)*nu/(h*h)-4/(tau*tau);
        M[2*m*i+1][n*m+i*m+1] = 2/(tau*tau);;
        M[2*m*i+1][n*m+(i+1)*m+m-1] = 0;
        M[2*m*i+1][n*m+(i+1)*m] = -(a+h*i)*(2*a*(nu-1)*nu+h*(-1+(-1+2*i*(nu-1))*nu))/(2*h*h);
        M[2*m*i+1][n*m+(i+1)*m+1] = 0;
        //
        M[2*m*i+2*(m-1)][(i-1)*m+m-2] = 0;
        M[2*m*i+2*(m-1)][(i-1)*m+m-1] = 2*a*a/(h*h) - (a/h) - i +4*a*i/h +2*i*i;
        M[2*m*i+2*(m-1)][(i-1)*m] = 0;
        M[2*m*i+2*(m-1)][i*m+m-2] = 1/(tau*tau) - (nu/(tau*tau));
        M[2*m*i+2*(m-1)][i*m+m-1] = -2-4*a*a/(h*h) - (8*a*i)/h - 4*i*i - 2/(tau*tau)+2*nu/(tau*tau);
        M[2*m*i+2*(m-1)][i*m] = 1/(tau*tau) - (nu/(tau*tau));
        M[2*m*i+2*(m-1)][(i+1)*m+m-2] = 0;
        M[2*m*i+2*(m-1)][(i+1)*m+m-1] = 2*a*a/(h*h) - (a/h) - i +4*a*i/h +2*i*i;
        M[2*m*i+2*(m-1)][(i+1)*m] = 0;
        M[2*m*i+2*(m-1)][n*m+(i-1)*m+m-2] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
        M[2*m*i+2*(m-1)][n*m+(i-1)*m+m-1] = 0;
        M[2*m*i+2*(m-1)][n*m+(i-1)*m] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
        M[2*m*i+2*(m-1)][n*m+i*m+m-2] = 3/(2*tau)-nu/(2*tau);
        M[2*m*i+2*(m-1)][n*m+i*m+m-1] = 0;
        M[2*m*i+2*(m-1)][n*m+i*m] = -3/(2*tau)+nu/(2*tau);
        M[2*m*i+2*(m-1)][n*m+(i+1)*m+m-2] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
        M[2*m*i+2*(m-1)][n*m+(i+1)*m+m-1] = 0;
        M[2*m*i+2*(m-1)][n*m+(i+1)*m] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
        //
        M[2*m*i+2*(m-1)+1][(i-1)*m+m-2] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
        M[2*m*i+2*(m-1)+1][(i-1)*m+m-1] = 0;
        M[2*m*i+2*(m-1)+1][(i-1)*m] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
        M[2*m*i+2*(m-1)+1][i*m+m-2] = -3/(2*tau)+nu/(2*tau);
        M[2*m*i+2*(m-1)+1][i*m+m-1] = 0;
        M[2*m*i+2*(m-1)+1][i*m] = 3/(2*tau)-nu/(2*tau);
        M[2*m*i+2*(m-1)+1][(i+1)*m+m-2] = a/(4*h*tau)+i/(4*tau)+a*nu/(4*h*tau)+i*nu/(4*tau);
        M[2*m*i+2*(m-1)+1][(i+1)*m+m-1] = 0;
        M[2*m*i+2*(m-1)+1][(i+1)*m] = -a/(4*h*tau)-i/(4*tau)-a*nu/(4*h*tau)-i*nu/(4*tau);
        M[2*m*i+2*(m-1)+1][n*m+(i-1)*m+m-2] = 0;
        M[2*m*i+2*(m-1)+1][n*m+(i-1)*m+m-1] = -(a+h*i)*(2*a*(nu-1)*nu+h*(-1+(-1+2*i*(nu-1))))/(2*h*h);
        M[2*m*i+2*(m-1)+1][n*m+(i-1)*m] = 0;
        M[2*m*i+2*(m-1)+1][n*m+i*m+m-2] = 2/(tau*tau);
        M[2*m*i+2*(m-1)+1][n*m+i*m+m-1] = -1+nu+2*(a+h*i)*(a+h*i)*(nu-1)*nu/(h*h)-4/(tau*tau);
        M[2*m*i+2*(m-1)+1][n*m+i*m] = 2/(tau*tau);;
        M[2*m*i+2*(m-1)+1][n*m+(i+1)*m+m-2] = 0;
        M[2*m*i+2*(m-1)+1][n*m+(i+1)*m+m-1] = -(a+h*i)*(2*a*(nu-1)*nu+h*(-1+(-1+2*i*(nu-1))*nu))/(2*h*h);
        M[2*m*i+2*(m-1)+1][n*m+(i+1)*m] = 0;
    }
    int jn = 0;
    for (int j = 0; 2*M_PI*j/m < M_PI/4; ++j) {
        M[2*m*(n-1)+2*j][(n-1)*m+j] = 1;
        B[2*m*(n-1)+2*j] = -U*std::cos(2*M_PI*j/m);
        M[2*m*(n-1)+2*j+1][n*m+(n-1)*m+j] = 1;
        B[2*m*(n-1)+2*j+1] = U*std::sin(2*M_PI*j/m);
        jn = j;
    }
    for (int j = jn+1; 2*M_PI*j/m < 3*M_PI/4; ++j) {
        M[2*m*(n-1)+2*j][(n-1)*m+j-2*m] = -E*b*tau/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j][(n-1)*m+j-m] = E*h*nu*tau/(h*b*tau-h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j][(n-1)*m+j] = E*b*tau/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j][n*m+(n-1)*m-m+j-1] = -E*h*nu/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j][n*m+(n-1)*m-m+j+1] = E*h*nu/(2*h*a*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][(n-1)*m+j-2*m] = -E*b*tau*nu/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][(n-1)*m+j-m] = E*h*tau/(h*b*tau-h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][(n-1)*m+j] = E*b*nu*tau/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][n*m+(n-1)*m-m+j-1] = -E*h/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][n*m+(n-1)*m-m+j+1] = E*h/(2*h*b*tau-2*h*b*nu*nu*tau);
        jn = j;
    }
    for (int j = jn+1; 2*M_PI*j/m < 5*M_PI/4; ++j){
        M[2*m*(n-1)+2*j][(n-1)*m+j] = 1;
        B[2*m*(n-1)+2*j] = U*std::cos(2*M_PI*j/m);
        M[2*m*(n-1)+2*j+1][n*m+(n-1)*m+j] = 1;
        B[2*m*(n-1)+2*j+1] = -U*std::sin(2*M_PI*j/m);
        jn = j;
    }
    for (int j = jn+1; 2*M_PI*j/m < 7*M_PI/4; ++j) {
        M[2*m*(n-1)+2*j][(n-1)*m+j-2*m] = -E*b*tau/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j][(n-1)*m+j-m] = E*h*nu*tau/(h*b*tau-h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j][(n-1)*m+j] = E*b*tau/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j][n*m+(n-1)*m-m+j-1] = -E*h*nu/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j][n*m+(n-1)*m-m+j+1] = E*h*nu/(2*h*a*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][(n-1)*m+j-2*m] = -E*b*tau*nu/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][(n-1)*m+j-m] = E*h*tau/(h*b*tau-h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][(n-1)*m+j] = E*b*nu*tau/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][n*m+(n-1)*m-m+j-1] = -E*h/(2*h*b*tau-2*h*b*nu*nu*tau);
        M[2*m*(n-1)+2*j+1][n*m+(n-1)*m-m+j+1] = E*h/(2*h*b*tau-2*h*b*nu*nu*tau);
        jn = j;
    }
    for (int j = jn+1; j<m; ++j) {
        M[2*m*(n-1)+2*j][(n-1)*m+j] = 1;
        B[2*m*(n-1)+2*j] = -U*std::cos(2*M_PI*j/m);
        M[2*m*(n-1)+2*j+1][n*m+(n-1)*m+j] = 1;
        B[2*m*(n-1)+2*j+1] = U*std::sin(2*M_PI*j/m);
        jn = j;
    }
    std::ofstream matrix_out("/home/vadik/matrix.txt");
    w_out(M,matrix_out);
    std::ofstream vector_out("/home/vadik/vector.txt");
    v_out(B,vector_out);
    return 0;
}
