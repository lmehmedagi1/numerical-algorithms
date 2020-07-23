//NA 2017/2018: Zadaća 1, Zadatak 1
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>
#include <iomanip>

class Vector {
    std::vector<double> v;
    
    public:
    
    explicit Vector(int n);
    Vector(std::initializer_list<double> l);
    int NElems() const;
    
    double &operator[](int i) {
        if (i<0 || i>=v.size()) throw std::range_error("Invalid index");
        return v.at(i); 
    }
    
    double operator[](int i) const {
        if (i<0 || i>=v.size()) throw std::range_error("Invalid index");
        return v.at(i);
    }
    
    double &operator()(int i) {
        if (i<1 || i>v.size()) throw std::range_error("Invalid index");
        return v.at(i-1);
    }
    
    double operator()(int i) const {
        if (i<1 || i>v.size()) throw std::range_error("Invalid index");
        return v.at(i-1);
    }
    
    double Norm() const;
    friend double VectorNorm(const Vector &v);
    double GetEpsilon() const;
    void Print(char separator = '\n', double eps = -1) const;
    friend void PrintVector(const Vector &v, char separator = '\n', double eps = -1);
    friend Vector operator +(const Vector &v1, const Vector &v2);
    
    Vector &operator +=(const Vector &v) {
        if (NElems()!=v.NElems()) throw std::domain_error("Incompatible formats");
        for (int i=0; i<NElems(); i++) Vector::v[i]+=v[i];
        return *this;
    }
    
    friend Vector operator -(const Vector &v1, const Vector &v2);
    
    Vector &operator -=(const Vector &v) {
        if (NElems()!=v.NElems()) throw std::domain_error("Incompatible formats");
        for (int i=0; i<NElems(); i++) Vector::v[i]-=v[i];
        return *this;
    }
    
    friend Vector operator *(double s, const Vector &v);
    friend Vector operator *(const Vector &v, double s);
    
    Vector &operator *=(double s) {
        for (double &x: v) x*=s;
        return *this;
    }
    
    friend double operator *(const Vector &v1, const Vector &v2);
    friend Vector operator /(const Vector &v, double s);
    
    Vector &operator /=(double s) {
        if (std::fabs(s)<std::numeric_limits<double>::epsilon()) throw std::domain_error("Division by zero");
        for (double &x: v) x/=s;
        return *this;
    }
};

Vector::Vector(int n) {
    if (n<=0) throw std::range_error("Bad dimension");
    for (int i=0; i<n; i++) v.push_back(0);
}

Vector::Vector(std::initializer_list<double> l) {
    if (l.size()==0) throw std::range_error("Bad dimension");
    for (double x: l) v.push_back(x);
}

int Vector::NElems() const {
    return v.size();
}

double Vector::Norm() const {
    double suma_kvadrata(0);
    for (int i=0; i<NElems(); i++) suma_kvadrata+=v.at(i)*v.at(i);
    return std::sqrt(suma_kvadrata);
}

double VectorNorm(const Vector &v) {
    return v.Norm();
}

double Vector::GetEpsilon() const {
    return 10*Norm()*std::numeric_limits<double>::epsilon();
}

void Vector::Print(char separator, double eps) const {
    if (eps<0) eps=GetEpsilon();
    
    for (int i=0; i<NElems(); i++) {
        if (std::fabs(v.at(i))<eps) std::cout << 0;
        else std::cout << v.at(i);
        if (i<NElems()-1 || separator=='\n') std::cout << separator;
    }
}
    
void PrintVector(const Vector &v, char separator, double eps) {
    v.Print(separator, eps);
}
    
Vector operator +(const Vector &v1, const Vector &v2) {
    if (v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector pomocni(v1.NElems());
    for (int i=0; i<v1.NElems(); i++) pomocni[i]=v1[i]+v2[i];
    return pomocni;
}
    
Vector operator -(const Vector &v1, const Vector &v2) {
    if (v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector pomocni(v1.NElems());
    for (int i=0; i<v1.NElems(); i++) pomocni[i]=v1[i]-v2[i];
    return pomocni;
}
    
Vector operator *(double s, const Vector &v) {
    Vector pomocni(v.NElems());
    for (int i=0; i<v.NElems(); i++) pomocni[i]=v[i]*s;
    return pomocni;
}

Vector operator *(const Vector &v, double s) {
    return s*v;
}
    
double operator *(const Vector &v1, const Vector &v2) {
    if (v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    double suma(0);
    for (int i=0; i<v1.NElems(); i++) suma+=v1[i]*v2[i];
    return suma;
}

Vector operator /(const Vector &v, double s) {
    if (std::fabs(s)<std::numeric_limits<double>::epsilon()) throw std::domain_error("Division by zero");
    Vector pomocni(v.NElems());
    for (int i=0; i<v.NElems(); i++) pomocni[i]=v[i]/s;
    return pomocni;
}


class Matrix {
    std::vector<std::vector<double>> mat;
    
    void ProvjeriJeLiNeispravnaLista(std::initializer_list<std::vector<double>> l) const {
        if (l.size()==0) throw std::range_error("Bad dimension");
        int br_kolona(0), brojac(0);
        for (std::vector<double> v: l) { br_kolona=v.size(); if (br_kolona==0) brojac++; }
        if (brojac) throw std::range_error("Bad dimension");
        for (std::vector<double> v: l) if (v.size()!=br_kolona) throw std::logic_error("Bad matrix");
    }
    
    public:
    
    Matrix(int m, int n) {
        if (m<=0 || n<=0) throw std::range_error("Bad dimension");
        mat.resize(m);
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++) {
                mat.at(i).push_back(0);
            }
        }
    }
    
    Matrix(const Vector &v) {
        mat.resize(v.NElems());
        for (int i=0; i<v.NElems(); i++) mat.at(i).push_back(v[i]);
    }
    
    Matrix(std::initializer_list<std::vector<double>> l) {
        ProvjeriJeLiNeispravnaLista(l);
        
        for (std::vector<double> v: l) mat.push_back(v); //?
    }
    
    int NRows() const { return mat.size(); }
    int NCols() const { return mat.at(0).size(); } //hmmmmmmmmm
    
    double *operator[](int i) {
        if (i<0 || i>=NRows()) throw std::range_error("Invalid index");
        return &mat.at(i).at(0);
    }
    
    const double *operator[](int i) const {
        if (i<0 || i>=NRows()) throw std::range_error("Invalid index");
        return &mat.at(i).at(0);
    }
    
    double &operator()(int i, int j) {
        if (i<1 || i>NRows() || j<1 || j>NCols()) throw std::range_error("Invalid index");
        return mat.at(i-1).at(j-1);
    }
    
    double operator()(int i, int j) const {
        if (i<1 || i>NRows() || j<1 || j>NCols()) throw std::range_error("Invalid index");
        return mat.at(i-1).at(j-1);
    }
    
    double Norm() const {
        double suma_kvadrata(0);
        for (int i=0; i<NRows(); i++) {
            for (int j=0; j<NCols(); j++) {
                suma_kvadrata+=mat.at(i).at(j)*mat.at(i).at(j);
            }
        }
        return std::sqrt(suma_kvadrata);
    }
    
    friend double MatrixNorm(const Matrix &m); 
    
    double GetEpsilon() const {
        return 10*Norm()*std::numeric_limits<double>::epsilon();
    }
    
    void Print(int width = 10, double eps = -1) const { // šta ako je width negativno
        if (eps<0) eps=GetEpsilon();
        for (int i=0; i<NRows(); i++) {
            for (int j=0; j<NCols(); j++) {
                if (std::fabs(mat.at(i).at(j))<eps) std::cout << std::setw(width) << 0;
                else std::cout << std::setw(width) << mat.at(i).at(j);
            }
            std::cout << std::endl;
        }
        
    }
    friend void PrintMatrix(const Matrix &m, int width = 10, double eps = -1);
    friend Matrix operator +(const Matrix &m1, const Matrix &m2);
    
    Matrix &operator +=(const Matrix &m) {
        if (NRows()!=m.NRows() || NCols()!=m.NCols()) throw std::domain_error("Incompatible formats");
        for (int i=0; i<NRows(); i++) {
            for (int j=0; j<NCols(); j++) {
                mat.at(i).at(j)+=m(i+1, j+1);
            }
        }
        return *this;
    }
    
    friend Matrix operator -(const Matrix &m1, const Matrix &m2);
    
    Matrix &operator -=(const Matrix &m) {
        if (NRows()!=m.NRows() || NCols()!=m.NCols()) throw std::domain_error("Incompatible formats");
        for (int i=0; i<NRows(); i++) {
            for (int j=0; j<NCols(); j++) {
                mat.at(i).at(j)-=m(i+1, j+1);
            }
        }
        return *this;
    }
    
    friend Matrix operator *(double s, const Matrix &m);
    friend Matrix operator *(const Matrix &m, double s);
    
    Matrix &operator *=(double s) {
        for (int i=0; i<NRows(); i++) {
            for (int j=0; j<NCols(); j++) {
                mat.at(i).at(j)*=s;
            }
        }
        return *this;
    }
    
    friend Matrix operator *(const Matrix &m1, const Matrix &m2);   
    
    Matrix &operator *=(const Matrix &m) {
        *this=*this*m;
        return *this;
    }
    
    friend Vector operator *(const Matrix &m, const Vector &v);
    
    friend Matrix Transpose(const Matrix &m);
    
    void Transpose() {
        if (NRows()==NCols()) {
            double pomocna;
            int x(1);
            for (int i=0; i<NRows(); i++) {
                for (int j=x; j<NCols(); j++) {
                    if (i!=j) {
                        pomocna=mat.at(i).at(j);
                        (*this)(i+1, j+1)=mat.at(j).at(i);
                        (*this)(j+1, i+1)=pomocna;
                    }
                }
                x++;
            }
            return;
        }
        
        Matrix pomocna(NCols(), NRows());
        for (int i=0; i<NCols(); i++) {
            for (int j=0; j<NRows(); j++) {
                pomocna(i+1, j+1)=mat.at(j).at(i);
            }
        }
        *this=pomocna;
    }

};

double MatrixNorm(const Matrix &m) {
    return m.Norm();
}

void PrintMatrix(const Matrix &m, int width, double eps) {
    m.Print(width, eps);
}

Matrix operator +(const Matrix &m1, const Matrix &m2) {
    if (m1.NRows()!=m2.NRows() || m1.NCols()!=m2.NCols()) throw std::domain_error("Incompatible formats");
    Matrix pomocna(m1.NRows(), m2.NCols());
    for (int i=0; i<m1.NRows(); i++) {
        for (int j=0; j<m2.NCols(); j++) {
            pomocna(i+1, j+1)=m1(i+1, j+1)+m2(i+1, j+1);
        }
    }
    return pomocna;
}

Matrix operator -(const Matrix &m1, const Matrix &m2) {
    if (m1.NRows()!=m2.NRows() || m1.NCols()!=m2.NCols()) throw std::domain_error("Incompatible formats");
    Matrix pomocna(m1.NRows(), m2.NCols());
    for (int i=0; i<m1.NRows(); i++) {
        for (int j=0; j<m2.NCols(); j++) {
            pomocna(i+1, j+1)=m1(i+1, j+1)-m2(i+1, j+1);
        }
    }
    return pomocna;
}

Matrix operator *(double s, const Matrix &m) {
    Matrix pomocna(m.NRows(), m.NCols());
    for (int i=0; i<m.NRows(); i++) {
        for (int j=0; j<m.NCols(); j++) {
            pomocna(i+1, j+1)=m(i+1, j+1)*s;
        }
    }
    return pomocna;
}

Matrix operator *(const Matrix &m, double s) { return s*m; }

Matrix operator *(const Matrix &m1, const Matrix &m2) {
    if (m1.NCols()!=m2.NRows()) throw std::domain_error("Incompatible formats");
    
    Matrix pomocna(m1.NRows(), m2.NCols());
    
    for (int i=0; i<m1.NRows(); i++) {
        for (int j=0; j<m2.NCols(); j++) {
            for(int k=0; k<m1.NCols(); k++) {
                pomocna(i+1, j+1)+=m1(i+1, k+1)*m2(k+1, j+1);
            }
        }
    }
    return pomocna;
}

Vector operator *(const Matrix &m, const Vector &v) {
    if (m.NCols()!=v.NElems()) throw std::domain_error("Incompatible formats");
    Vector pomocni(m.NRows());
    
    for (int i=0; i<m.NRows(); i++) {
        for(int k=0; k<m.NCols(); k++) {
            pomocni(i+1)+=m(i+1, k+1)*v(k+1);
        }
    }
    return pomocni;
}

Matrix Transpose(const Matrix &m) {
    Matrix pomocna(m.NCols(), m.NRows());
    for (int i=0; i<m.NCols(); i++) {
        for (int j=0; j<m.NRows(); j++) {
            pomocna(i+1, j+1)=m(j+1, i+1);
        }
    }
    return pomocna;
}



int main ()
{
    //TESTIRANJE KLASE VECTOR
    
    //1. Konstruktor sa cjelobrojnim parametrom
    Vector v1(5);
    Vector v2(5);
    
    try {
        Vector v3(-1);
    }
    catch(std::range_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    
    
    //2. Sekvencijski konstruktor
    Vector v4({4, 4, 4, 4});
    try {
        Vector v5({});
    }
    catch(std::range_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    
    //3. NElems, operator [] i operator ()
    for (int i=0; i<v4.NElems(); i++) {
        v4[i]++;
    }
    for (int i=0; i<v1.NElems(); i++) {
        v1(i+1)=i*2;
    }
    try {
        std::cout << v1[7];
    }
    catch(std::range_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    
    std::cout << std::endl;
    
    //4. Print i PrintVector
    const Vector v(v4);
    v4.Print(',', 1); std::cout << std::endl;
    PrintVector(v4, ' ', -8);
    
    //5. Norm, VectorNorm i GetEpsilon
    if (std::fabs(v4.Norm()-VectorNorm(v4))>v4.GetEpsilon()) std::cout << "Greška!" << std::endl;
    
    //6. Operatori +, -, += i -=
    Vector v5(v1+v2);
    v2=v1-v2;
    try {
        v2+=v4;
    }
    catch(std::domain_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    v2-=v1;
    
    
    //7. Operatori *, *=, / i /=
    v5=v1*6;
    v5/=6;
    try {
        Vector v7(v4/0);
    }
    catch(std::domain_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    try {
        Vector v8(v1*v4);
    }
    catch(std::domain_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    
    v5*=3;
    
    //TESTIRANJE KLASE MATRIX
    
    //1. Konstruktor sa dva cjelobrojna parametra
    Matrix m1(2, 3);
    Matrix m2(3, 3);
    Matrix m3(2, 3);
    try {
        Matrix m4(0, 1);
    }
    catch(std::range_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    
    //2. Konstruktor koji prima vektor
    Matrix m4(v1);
    
    //3. Sekvencijski konstruktor
    const Matrix m5({{1, 2, 3}, {3, 3, 3}, {5, 6, 7}, {2, 1, 2}});
    try {
        Matrix m6({});
    }
    catch(std::range_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    try {
        Matrix m6({{1, 2}, {3, 2}, {}});
    }
    catch(std::range_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    try {
        Matrix m6({{1, 1, 1, 1}, {1, 1, 1}, {1, 1, 1, 1}});
    }
    catch(std::logic_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    
    //4. NRows, NCols, operator () i operator []
    for (int i=0; i<m1.NRows(); i++) {
        for (int j=0; j<m1.NCols(); j++) {
            m1(i+1, j+1)=2*i+j;
        }
    }
    std::cout << m5[2][1] << std::endl;
    try {
        m2(8,1)++;
    }
    catch(std::range_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    
    //5. Print i PrintMatrix
    m4.Print(4 , 1);
    PrintMatrix(m4, 4, -7);
    
    //6. Norm, MatrixNorm i GetEpsilon
    if (std::fabs(m1.Norm()-MatrixNorm(m1))>m1.GetEpsilon()) std::cout << "Greška";
    
    //7. Operatori +, -, += i -=
    Matrix m7(m1+m3);
    Matrix m8(m1-m3);
    m7+=m8;
    m8-=m7;
    try {
        PrintMatrix(m1-m2);
    }
    catch(std::domain_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    try {
        std::cout << MatrixNorm(m1+=m2);
    }
    catch(std::domain_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    
    //8. Operatori * i *=
    PrintMatrix(m1*4);
    PrintMatrix(3*m4);
    m7*=2;
    PrintMatrix(m1*m2);
    try {
        PrintMatrix(m2*v4);
    }
    catch(std::domain_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    try {
        std::cout << MatrixNorm(m2*m8);
    }
    catch(std::domain_error &iz) {
        std::cout << iz.what() << std::endl;
    }
    
    //9. Transpose
    std::cout << "Prije transponovanja: " << std::endl;
    m1.Print(4);
    std::cout << std::endl << "Nakon transponovanja: " << std::endl;
    m1.Transpose();
    m1.Print(4);
    
    
	return 0;
}