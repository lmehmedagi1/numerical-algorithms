//NA 2018/2019: Zadaća 2, Zadatak 1
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
    
    void Chop(double eps=-1) {
        if (eps<0) eps=GetEpsilon();
        for (int i=0; i<v.size(); i++) {
            if (std::fabs(v.at(i))<eps) v.at(i)=0;
        }
    }
    
    bool EqualTo(const Vector &vek, double eps=-1) const {
        if (eps<0) eps=GetEpsilon();
        if (v.size()!=vek.v.size()) return false;
        for (int i=0; i<v.size(); i++) if (std::fabs(v.at(i)-vek.v.at(i))>eps) return false;
        return true;
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
    
    void ZamijeniKolone(int i, int j) {
        double temp;
        for (int k=0; k<NRows(); k++) {
            temp=mat[k][i];
            mat[k][i]=mat[k][j];
            mat[k][j]=temp;
        }
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
    
    
    void Chop(double eps = -1) {
        if (eps<0) eps=GetEpsilon();
        for (int i=0; i<NRows(); i++) {
            for (int j=0; j<NCols(); j++) {
                if (std::fabs(mat.at(i).at(j))<eps) mat.at(i).at(j)=0;
            }
        }
        
    }
    
    bool EqualTo(const Matrix &m, double eps = -1) const {
        if (eps<0) eps=GetEpsilon();
        if (NCols()!=m.NCols() || NRows()!=m.NRows()) return false;
        for (int i=0; i<NRows(); i++) {
            for (int j=0; j<NCols(); j++) {
                if (std::fabs(mat.at(i).at(j)-m.mat.at(i).at(j))>eps) return false;
            }
        }
        return true;
    }
    
    friend Matrix LeftDiv(Matrix m1, Matrix m2);
    
    friend Vector LeftDiv(Matrix m, Vector v);
    
    friend Matrix operator /(const Matrix &m, double s);
    
    Matrix &operator /=(double s) {
        if (std::fabs(s)<GetEpsilon()) throw std::domain_error("Division by zero");
        for (int i=0; i<NRows(); i++) {
            for (int j=0; j<NCols(); j++) {
                mat.at(i).at(j)/=s;
            }
        }
        return *this;
    }
    
    friend Matrix operator /(Matrix m1, Matrix m2);
    
    
    Matrix &operator /=(Matrix m) {
        if (m.NCols()!=m.NRows()) throw std::domain_error("Divisor matrix is not square");
        if (std::fabs(m.Det())<GetEpsilon()) throw std::domain_error("Divisor matrix is singular");
        if (NCols()!=m.NCols()) throw std::domain_error("Incompatible formats");
        
        double mi(0);
        int p;
        
        Matrix x(NRows(), m.NCols());
        
        for (int k=0; k<m.NRows(); k++) {
            p=k;
            for (int i=k+1; i<m.NRows(); i++) {
                if (std::fabs(m[k][i])>std::fabs(m[k][p])) p=i;
            }
            if (std::fabs(m[k][p])<m.GetEpsilon()) throw std::domain_error("Divisor matrix is singular");
            if (p!=k) {
                m.ZamijeniKolone(p, k);
                ZamijeniKolone(p, k);
            }
            
            for (int i=k+1; i<m.NRows(); i++) {
                mi=m[k][i]/m[k][k];
                for (int j=k+1; j<m.NRows(); j++) {
                    m[j][i]=m[j][i]-mi*m[j][k];
                }
                for (int j=0; j<NRows(); j++) {
                    (*this)[j][i]=(*this)[j][i]-mi*(*this)[j][k];
                }
            }
        }
        
        
        double s(0);
        for (int k=0; k<NRows(); k++) {
            for (int i=m.NRows()-1; i>=0; i--) {
                s=(*this)[k][i];
                for (int j=i+1; j<m.NRows(); j++) {
                    s=s-m[j][i]*x[k][j];
                }
                x[k][i]=s/m[i][i];
            }
        }
        *this=x;
        return *this; 
    }
    
    
    friend double Det(Matrix m);
    
    
    double Det() const {
        if (NCols()!=NRows()) throw std::domain_error("Matrix is not square");
        Matrix m(*this);
        if (m.NCols()!=m.NRows()) throw std::domain_error("Matrix is not square");
        double d(1), mi(0);
        int p;
        for (int k=0; k<m.NRows(); k++) {
            p=k;
            for (int i=k+1; i<m.NRows(); i++) {
                if (std::fabs(m[i][k])>std::fabs(m[p][k])) p=i;
            }
            if (std::fabs(m[p][k])<m.GetEpsilon()) return 0;
            if (p!=k) {
                (m.mat.at(k)).swap(m.mat.at(p));
                d=-d;
            }
            d=d*m[k][k];
            for (int i=k+1; i<m.NRows(); i++) {
                mi=m[i][k]/m[k][k];
                for (int j=k+1; j<m.NRows(); j++) {
                    m[i][j]=m[i][j]-mi*m[k][j];
                }
            }
        }
        return d;
    }
    
    void ZamijeniRedove(int i, int j) {
        mat.at(i).swap(mat.at(j));
    }
    
    
    void Invert() {
        if (NRows()!=NCols()) throw std::domain_error("Matrix is not square");
        if (std::fabs(Det())<GetEpsilon()) throw std::domain_error("Matrix is singular");
        
        std::vector<int> w; w.resize(NRows());
        double mi;
        int p;
        
        for (int k=0; k<NRows(); k++) {
            p=k;
            for (int i=k+1; i<NRows(); i++) {
                if (std::fabs((*this)[i][k])>std::fabs((*this)[p][k])) p=i;
            }
            if (std::fabs((*this)[p][k])<GetEpsilon()) throw std::domain_error("Matrix is singular");
            if (p!=k) {
                mat.at(k).swap(mat.at(p));
            }
            w.at(k)=p;
            mi=(*this)[k][k];
            (*this)[k][k]=1;
            for (int j=0; j<NRows(); j++) {
                (*this)[k][j]=(*this)[k][j]/mi;
            }
            for (int i=0; i<NRows(); i++) {
                if (i!=k) {
                    mi=(*this)[i][k];
                    (*this)[i][k]=0;
                    for (int j=0; j<NRows(); j++) {
                        (*this)[i][j]=(*this)[i][j]-mi*(*this)[k][j];
                    }
                }
            }
        }
        
        for(int j=NRows()-1; j>=0; j--) {
            p=w.at(j);
            if (p!=j) {
                for (int i=0; i<NRows(); i++) {
                    ZamijeniKolone(p, j);
                }
            }
        }
    }
    
    
    friend Matrix Inverse(Matrix m);
    
    void ReduceToRREF() {
        
        int k(-1), l(-1), p(0);
        double mi, v;
        std::vector<bool> w; w.resize(NCols());
        
        for (int j=0; j<NCols(); j++) w.at(j)=false;
        
        while (k<NRows() && l<NCols()) {
            l++;
            k++;
            v=0;
            while (v<GetEpsilon() && l<NCols()) {
                p=k;
                for (int i=k; i<NRows(); i++) {
                    if (std::fabs((*this)[i][l])>v) {
                        v=std::fabs((*this)[i][l]);
                        p=i;
                    }
                }
                if (v<GetEpsilon()) l++;
            }
            if (l<NCols()) {
                w.at(l)=true;
                if (p!=k) {
                    mat.at(k).swap(mat.at(p));
                }
                mi=(*this)[k][l];
                for (int j=l; j<NCols(); j++) {
                    (*this)[k][j]=(*this)[k][j]/mi;
                }
                for (int i=0; i<NRows(); i++) {
                    if (i!=k) {
                        mi=(*this)[i][l];
                        for (int j=l; j<NCols(); j++) {
                            (*this)[i][j]=(*this)[i][j]-mi*(*this)[k][j];
                        }
                    }
                    
                }
            }
        }
    }
    
    friend Matrix RREF(Matrix m);
    
    int Rank() const {
        Matrix pomocna(*this);
        pomocna.ReduceToRREF();
        int rang(0);
        for (int i=0; i<pomocna.NRows(); i++) {
            for (int j=0; j<pomocna.NCols(); j++) {
                if (std::fabs(pomocna[i][j])>GetEpsilon()) { rang++; break; }
            }
        }
        return rang;
    }
    
    friend int Rank(Matrix m);

};

int Rank(Matrix m) {
    return m.Rank();
}

Matrix RREF(Matrix m) {
    m.ReduceToRREF();
    return m;
}

Matrix Inverse(Matrix m) {
    m.Invert();
    return m;
}

Matrix operator /(const Matrix &m, double s) {
    Matrix pomocna(m);
    return pomocna/=s;
}

Matrix operator /(Matrix m1, Matrix m2) {
    return m1/=m2;
}

Matrix LeftDiv(Matrix m1, Matrix m2) {
    if (m1.NCols()!=m1.NRows()) throw std::domain_error("Divisor matrix is not square");
    if (m1.NRows()!=m2.NRows()) throw std::domain_error("Incompatible formats");
    if (std::fabs(m1.Det())<m1.GetEpsilon()) throw std::domain_error("Divisor matrix is singular");
    
    double mi(0);
    int p;
    
    Matrix x(m1.NRows(), m2.NCols());
    
    for (int k=0; k<m1.NRows(); k++) {
        p=k;
        for (int i=k+1; i<m1.NRows(); i++) {
            if (std::fabs(m1[i][k])>std::fabs(m1[p][k])) p=i;
        }
        if (std::fabs(m1[p][k])<m1.GetEpsilon()) throw std::domain_error("Divisor matrix is singular");
        if (p!=k) {
            m1.mat.at(k).swap(m1.mat.at(p));
            m2.mat.at(k).swap(m2.mat.at(p));
        }
        
        for (int i=k+1; i<m1.NRows(); i++) {
            mi=m1[i][k]/m1[k][k];
            for (int j=k+1; j<m1.NRows(); j++) {
                m1[i][j]=m1[i][j]-mi*m1[k][j];
            }
            for (int j=0; j<m2.NCols(); j++) {
                m2[i][j]=m2[i][j]-mi*m2[k][j];
            }
        }
    }
    
    
    double s(0);
    for (int k=0; k<m2.NCols(); k++) {
        for (int i=m1.NRows()-1; i>=0; i--) {
            s=m2[i][k];
            for (int j=i+1; j<m1.NRows(); j++) {
                s=s-m1[i][j]*x[j][k];
            }
            x[i][k]=s/m1[i][i];
        }
    }
    return x;
}

Vector LeftDiv(Matrix m, Vector v) {
    Matrix pomocna(LeftDiv(m, Matrix(v)));
    Vector vek(pomocna.NRows());
    for (int i=0; i<pomocna.NRows(); i++) vek[i]=pomocna.mat.at(i).at(0);
    return vek;
}

double MatrixNorm(const Matrix &m) {
    return m.Norm();
}

void PrintMatrix(const Matrix &m, int width, double eps) {
    m.Print(width, eps);
}

double Det(Matrix m) {
        return m.Det();
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



class LUDecomposer {
    
    Matrix a;
    Vector w;
    
    public:
    LUDecomposer(Matrix m) : a(m), w(m.NRows()) {
        if (m.NRows()!=m.NCols()) throw std::domain_error("Matrix is not square");
        if (std::fabs(m.Det())<m.GetEpsilon()) throw std::domain_error("Matrix is singular"); 
        
        int p;
        double s, mi;
        
        for (int j=0; j<a.NRows(); j++) {
            for (int i=0; i<=j; i++) {
                s=a[i][j];
                for (int k=0; k<=i-1; k++) {
                    s=s-a[i][k]*a[k][j];
                }
                a[i][j]=s;
            }
            p=j;
            for (int i=j+1; i<a.NRows(); i++) {
                s=a[i][j];
                for (int k=0; k<=j-1; k++) {
                    s=s-a[i][k]*a[k][j];
                }
                a[i][j]=s;
                if (std::fabs(s)>std::fabs(a[p][j])) {
                    p=i;
                }
            }
            if (std::fabs(a[p][j])<a.GetEpsilon()) throw std::domain_error("Matrix is singular");
            if (p!=j) a.ZamijeniRedove(j, p);
            w[j]=p;
            mi=a[j][j];
            for (int i=j+1; i<a.NRows(); i++) a[i][j]=a[i][j]/mi;
        }
    }
    
    void Solve(const Vector &b, Vector &x) const {
        if(b.NElems()!=a.NRows() || x.NElems()!=a.NRows()) throw std::domain_error("Incompatible formats");
        
        int p;
        double s;
        
        x=b;
        
        for (int i=0; i<a.NRows(); i++) {
            p=w[i];
            s=x[p];
            x[p]=x[i];
            for (int j=0; j<=i-1; j++) {
                s=s-a[i][j]*x[j];
            }
            x[i]=s;
        }
        
        for(int i=b.NElems()-1; i>=0; i--) {
            s=x[i];
            for(int j=i+1; j<b.NElems(); j++) {
                s=s-a[i][j]*x[j];
            }
            x[i]=s/a[i][i];
        }
    }
    
    Vector Solve(Vector b) const {
        Vector x(b.NElems());
        Solve(b, x);
        return x;
    }
    
    void Solve(Matrix &b, Matrix &x) const {
        if(b.NCols()!=a.NRows() || b.NRows()!=a.NCols() || x.NRows()!=a.NCols() || x.NCols()!=a.NRows()) throw std::domain_error("Incompatible formats");
        x=b;
        for(int k=0; k<b.NRows(); k++) {
            for(int i=0; i<a.NCols(); i++) {
                int p=w[i];
                double s=x[p][k];
                x[p][k]=x[i][k];
                for(int j=0; j<i; j++) {
                    s=s-a[i][j]*x[j][k];
                }
                x[i][k]=s;
            }
            for(int i(a.NCols()-1);i>=0;i--) {
                double s=x[i][k];
                for(int j(i+1);j<a.NCols();j++) {
                    s-=a[i][j]*x[j][k];
                }
                x[i][k]=s/a[i][i];
            }   
        }
    }
    
    Matrix Solve(Matrix b) const {
        Matrix x(b.NRows(), b.NCols());
        Solve(b, x);
        return x;
    }   
    
    Matrix GetCompactLU() const {
        return a;
    }
    
    Matrix GetL() const {
        Matrix l(a.NRows(), a.NCols());
        for (int i=0; i<a.NRows(); i++) {
            for(int j=0; j<a.NCols(); j++) {
                if (i==j) l[i][j]=1;
                else if (i>j) l[i][j]=a[i][j];
                else l[i][j]=0;
            } 
        }
        return l;
    }
    
    Matrix GetU() const {
        Matrix u(a.NRows(), a.NCols());
        for (int i=0; i<a.NRows(); i++) {
            for(int j=0; j<a.NCols(); j++) {
                if (i<=j) u[i][j]=a[i][j];
                else u[i][j]=0;
            }
        }
        return u;
    }
    
    Vector GetPermuation() const {
        return w;
    }
};


class QRDecomposer {
    Matrix a;
    std::vector<double> d;
    
    public:
    QRDecomposer(Matrix m) : a(m) {
        if (m.NRows()<m.NCols()) throw std::domain_error("Invalid matrix format");
        
        d.resize(a.NCols());
        double mi(0), s(0);
        
        for (int k=0; k<a.NCols(); k++) {
           s=0;
           for (int i=k; i<a.NRows(); i++) {
               s=s+a[i][k]*a[i][k];
           }
           s=std::sqrt(s);
           mi=std::sqrt(s*(s+std::fabs(a[k][k])));
           if(std::fabs(mi)<a.GetEpsilon()) throw std::domain_error("Matrix is singular");
           if (a[k][k]<0) s=-s;
           a[k][k]=(a[k][k]+s)/mi;
           for (int i=k+1; i<a.NRows(); i++) {
               a[i][k]=a[i][k]/mi;
           }
           d.at(k)=-s;
           for (int j=k+1; j<a.NCols(); j++) {
               s=0;
               for (int i=k; i<a.NRows(); i++) {
                   s=s+a[i][k]*a[i][j];
               }
               for (int i=k; i<a.NRows(); i++) {
                   a[i][j]=a[i][j]-s*a[i][k];
               }
           }
        }
    }
    
    void Solve(const Vector &b, Vector &x) const {
        if (a.NRows()!=a.NCols()) throw std::domain_error("Matrix is not square");
        if(b.NElems()!=a.NRows() || b.NElems()!=x.NElems()) throw std::domain_error("Incompatible formats");
        Vector b2(MulQTWith(b));
        
        for(int i=a.NCols()-1; i>=0; i--) {
            double s=b2[i];
            for(int j=i+1; j<a.NCols(); j++) {
                
                if (j>i) s=s-a[i][j]*x[j];
                else if (i==j) s=s-d[i]*x[j];
            }
            x[i]=s/d[i];
        }
    }
    
    Vector Solve(Vector b) const {
        Vector v(b.NElems());
        Solve(b,v);
        return v;
    }
    
    void Solve(Matrix &b, Matrix &x) const {
        if(a.NRows()!=a.NCols()) throw std::domain_error("Matrix is not square");
        if(b.NCols()!=a.NRows() || b.NRows()!=a.NCols() || x.NRows()!=a.NCols() || x.NCols()!=a.NRows()) throw std::domain_error("Incompatible formats");
        
        for(int i=0; i<b.NCols(); i++) {
            Vector v1(b.NRows()), v2(b.NRows());
            for(int j=0; j<b.NRows(); j++) {
                v1[j]=b[j][i];
            }
            Solve(v1,v2);
            for(int j=0; j<b.NRows(); j++) {
                x[j][i]=v2[j];
            }
        }
    }
        
    Matrix Solve(Matrix b) const {
        Matrix m(b.NRows(), b.NCols());
        Solve(b,m);
        return m;
    }
    
    Vector MulQWith(Vector v) const {
        if(a.NRows()!=v.NElems()) throw std::domain_error("Incompatible formats"); //hm
        
        double s;
        for (int k=a.NCols()-1; k>=0; k--) {
            s=0;
            for (int i=k; i<a.NRows(); i++) {
                s=s+a[i][k]*v[i];
            }
            for (int i=k; i<a.NRows(); i++) {
                v[i]=v[i]-s*a[i][k];
            }
        }
        return v;
    }
    
    Matrix MulQWith(Matrix m) const {
        if(a.NRows()!=m.NRows()) throw std::domain_error("Incompatible formats");
        
        for(int i=0; i<m.NCols(); i++) {
            Vector v(a.NRows());
            for(int j=0; j<v.NElems(); j++) v[j]=m[j][i];
            v=MulQWith(v);
            for(int j=0; j<v.NElems(); j++) m[j][i]=v[j];
        }
        return m;
    }
    
    Vector MulQTWith(Vector v) const {
        if(a.NRows()!=v.NElems()) throw std::domain_error("Incompatible formats"); 
        
        double s;
        for (int k=0; k<a.NCols(); k++) {
            s=0;
            for (int i=k; i<a.NRows(); i++) {
                s=s+a[i][k]*v[i];
            }
            for (int i=k; i<a.NRows(); i++) {
                v[i]=v[i]-s*a[i][k];
            }
        }
        return v;
    }
    
    Matrix MulQTWith(Matrix m) const {
        if(a.NRows()!=m.NRows()) throw std::domain_error("Incompatible formats");
        
        for(int i=0; i<m.NCols(); i++) {
            Vector v(a.NRows());
            for(int j=0; j<v.NElems(); j++) v[j]=m[j][i];
            v=MulQTWith(v);
            for(int j=0; j<v.NElems(); j++) m[j][i]=v[j];
        }
        return m;
    }
    
    Matrix GetQ() const {
        Matrix q(a.NRows(), a.NRows());
        for(int j=0; j<a.NRows(); j++) {
            for(int i=0; i<a.NRows(); i++) q[i][j]=0;
            q[j][j]=1;
            for(int k=a.NCols()-1; k>=0; k--) {
                double s=0;
                for(int i=k; i<a.NRows(); i++) {
                    s=s+a[i][k]*q[i][j];
                }
                for(int i=k; i<a.NRows(); i++) {
                    q[i][j]=q[i][j]-s*a[i][k];
                }
            }
        }
        return q;
    }
    
    Matrix GetR() const {
        Matrix r(a.NRows(),a.NCols());
        for(int i=0; i<r.NRows(); i++) {
            for(int j=0; j<r.NCols(); j++) { 
                if (j>i) r[i][j]=a[i][j];
                else if (j==i) r[i][j]=d[i];
                else r[i][j]=0;
            }
        }
        return r;
    }

};


//TESTOVI ZA VEKTOR:

void VectorChopTest() {
    Vector v({1, -2, 4, -5, 5, -4, 2, 1});
    v.Chop(3);
    if (std::fabs(v[1])>0 || std::fabs(v[3])<3) std::cout << "Test Vector Chop: Greska" << std::endl;
    else std::cout << "Test Vector Chop: OK" << std::endl;
}

void VectorEqualToTest() {
    Vector v1({1, -2, 4, -5, 5, -4, 2, 1});
    Vector v2({1, -2.5, 3.9, -5, 5, -4.1, 2, 1});
    Vector v3({1, -2, 4, 5, 5, -4, 2, 1});
    Vector v4({1, -2, 4, -5, 5, -4, 2});
    if (v1.EqualTo(v2, 0.5) && !v1.EqualTo(v3) && !v1.EqualTo(v4)) std::cout << "Test Vector EqualTo: OK" << std::endl;
    else std::cout << "Test Vector EqualTo: Greška" << std::endl;
}

//TESTOVI ZA MATRICU:

void MatrixChopTest() {
    Matrix m({{1, 2}, {2, 2}, {3, -3}});
    m.Chop(1.1);
    if (std::fabs(m[0][0])>0 || std::fabs(m[2][1])<3) std::cout << "Test Matrix Chop: Greska" << std::endl;
    else std::cout << "Test Matrix Chop: OK" << std::endl;
}

void MatrixEqualToTesst() {
    Matrix m1({{1, 2}, {2, 2}, {3, -3}});
    Matrix m2({{1, 2.1}, {1.9, 2}, {3, -2.8}});
    Matrix m3({{1, 2}, {2, 2}, {3, 3}});
    Matrix m4({{1, 2, 1}, {2, 2, 1}, {3, -3, 1}});
    if (m1.EqualTo(m2, 0.5) && !m1.EqualTo(m3) && !m3.EqualTo(m4)) std::cout << "Test Matrix EqualTo: OK" << std::endl;
    else std::cout << "Test Matrix EqualTo: Greška" << std::endl;
}

void LeftDivTest() {
    int brojac(0);
    Matrix m1({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Matrix m2({{1}, {2}, {3}});
    Matrix m3(LeftDiv(m1, m2));
    
    if (m3.NRows()==3 && m3.NCols()==1) brojac++;
    
    //Incompatible formats:
    try {
        Matrix m4({{1, 1, 1}, {2, 2, 2}});
        LeftDiv(m1,m4);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    //Divisor matrix is not square:
    try {
        Matrix m4({{1, 1, 1}, {2, 2, 2}});
        Matrix m5({{1, 1}, {2, 2}});
        LeftDiv(m4,m5);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    //Divisor matrix is singular:
    try {
        Matrix m4({{1, 0, 0}, {-2, 0, 0},{4, 6, 1}});
        Matrix m5({{1, 1}, {2, 2}, {3, 3}});
        LeftDiv(m4,m5);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    //Dijelenje vektorom:
    Matrix m({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Vector v({5, 5, 5});
    Vector rez(LeftDiv(m, v));
    
    if (v.NElems()==3) brojac++;
    
    if (brojac==5) std::cout << "Test LeftDiv: OK" << std::endl;
    else std::cout << "Test LeftDiv: Greška" << std::endl;
}

void RightDivTest() {
    int brojac(0);
    Matrix m{{3, 6, 9}, {12, 12, 12}};
    Matrix m1(m/3);
    Matrix rez({{1, 2, 3},{4, 4, 4}});
    if (m1.EqualTo(rez)) brojac++;
    Matrix M(m);
    M/=3;
    if (M.EqualTo(rez)) brojac++;
    
    //Division by zero:
    try {
        m1/0;
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    
    Matrix m2({{1, 2, 1},{1, 1, 3},{1, 1, 1}});
    Matrix m3(m/m2);
    if (m3.NRows()==2 && m3.NCols()==3) brojac++;
    
    
    //Incompatible formats:
    Matrix m4({{4, 5},{5, 4},{4, 5}});
    try {
        m4/m2;
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    //Divisor matrix is not square:
    try {
        m2/m4;
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    //Divisor matrix is singular:
    try {
        Matrix m4({{1, 0, 0}, {-2, 0, 0},{4, 6, 1}});
        Matrix m5({{1, 1}, {2, 2}, {3, 3}});
        LeftDiv(m5,m4);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    if (brojac==7) std::cout << "Test RightDiv: OK" << std::endl;
    else std::cout << "Test RightDiv: Greška" << std::endl;
}

void DetTest() {
    int brojac(0);
    Matrix m({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    if (std::fabs(Det(m)+37)<m.GetEpsilon()) brojac++;
    Matrix sing({{1, 0, 0}, {-2, 0, 0},{4, 6, 1}});
    if (std::fabs(sing.Det())<m.GetEpsilon()) brojac++;
    if (brojac==2) std::cout << "Test Det: OK" << std::endl;
    else std::cout << "Test Det: Greška" << std::endl;
}

void InvertTest() {
    int brojac(0);
    Matrix m({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Matrix rez({{-1.6486486, 1.4864865, -0.1081081}, {0.8378378, -0.3783784, -0.027027 }, {0.3243243, -0.2432432, 0.0540541}});
    if (Inverse(m).EqualTo(rez, 0.5)) brojac++;
    m.Invert();
    if (m.EqualTo(rez, 0.5)) brojac++;
    
    //Singular matrix:
    try {
        Matrix m4({{1, 0, 0}, {-2, 0, 0},{4, 6, 1}});
        Inverse(m4);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    //Matrix is not square:
    try {
        Matrix m4({{1, 0}, {-2, 0},{6, 1}});
        Inverse(m4);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    
    if (brojac==4) std::cout << "Test invert: OK" << std::endl;
    else std::cout << "Test invert: Greška" << std::endl;
}

void RREFTest() {
    int brojac(0);
    Matrix A{{0,3,2},{4,6,1},{3,1,7}};
    Matrix rez({{1,0,0}, {0,1,0}, {0,0,1}});
    if (RREF(A).EqualTo(rez)) brojac++;
    A.ReduceToRREF();
    if (A.EqualTo(rez)) brojac++;
    if (Rank(A)==3) brojac++;
    if (A.Rank()==3) brojac++;
    
    if (brojac==4) std::cout << "Test RREF: OK" << std::endl;
    else std::cout << "Test RREF: Greška" << std::endl;
}

//TESTOVI ZA LU:

void LUDecomposerTest() {
    int brojac(0);
    Matrix m({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Matrix rezL({{1, 0, 0}, {0.6666667, 1, 0}, {0.3333333, 0.75, 1}});
    Matrix rezU({{3, -3, 23}, {0, 4, -10.333333}, {0, 0, 3.0833333}});
    LUDecomposer lu(m);
    if (lu.GetL().EqualTo(rezL, 0.5)) brojac++;
    if (lu.GetU().EqualTo(rezU, 0.01)) brojac++;
    
    
    Vector rowswap({2, 1, 2});
    if (lu.GetPermuation().EqualTo(rowswap, 0.01)) brojac++;
    
    
    Matrix LiU({{3, -3, 23}, {0.666667, 4, -10.3333}, {0.333333, 0.75, 3.08333}});
    if (lu.GetCompactLU().EqualTo(LiU, 0.01)) brojac++;
    
    
    //Singular matrix:
    try {
        Matrix m4({{1, 0, 0}, {-2, 0, 0},{4, 6, 1}});
        LUDecomposer lu1(m4);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    
    //Matrix is not square:
    try {
        Matrix m4({{1, 0}, {-2, 0},{6, 1}});
        LUDecomposer lu1(m4);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    
    if (brojac==6) std::cout << "Test LUDecomposer: OK" << std::endl;
    else std::cout << "Test LUDecomposer: Greška" << std::endl;
}

void SolveLUTest() {
    int brojac(0);
    Matrix A({{2, 4}, {3, -3}});
    LUDecomposer lu(A);
    Vector b({-2, 5});
    Vector x(2);
    lu.Solve(b, x);
    if (x.EqualTo({0.777778, -0.888889}, 0.0001)) brojac++;
    if (lu.Solve(b).EqualTo({0.777778, -0.888889}, 0.0001)) brojac++;
    
    //Incompatible formats:
    Vector b2({1, 2, 3});
    try {
        lu.Solve(b2, x);
    }
    catch(std::domain_error) {brojac++;}
    catch(...) {}
    
    //LUSolve sa matricama:
    Matrix m1({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Matrix b1({{4, 4, 4}, {8, 7, 1}, {2, 2, 3}});
    Matrix x1(3,3);
    LUDecomposer lu1(m1);
    lu1.Solve(b1, x1);
    
    Matrix xrez({{5.08108, 3.59459, -5.43243}, {0.27027, 0.648649, 2.89189}, {-0.540541, -0.297297, 1.21622}});
    if (x1.EqualTo(xrez, 0.0001)) brojac++;
    if (lu1.Solve(b1).EqualTo(xrez, 0.0001)) brojac++;
    
    if (brojac==5) std::cout << "Test Solve LU: OK" << std::endl;
    else std::cout << "Test Solve LU: Greška" << std::endl;
}


//TESTOVI ZA QR:

void QRDecomposerTest() {
    int brojac(0);
    Matrix m({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    QRDecomposer qr(m);
    
    Matrix rrez({{-3.7416574, 0.8017837, -21.915422}, {0, -4.0443965, 8.7599155}, {0, 0, 2.4450288}});
    Matrix qrez({{-0.2672612, -0.5474947, 0.7929823}, {-0.5345225, -0.6004781, -0.5947367}, {-0.8017837, 0.582817, 0.1321637}});
    
    if (qr.GetR().EqualTo(rrez, 0.001)) brojac++;
    if (qr.GetQ().EqualTo(qrez, 0.001)) brojac++;
    
    try {
        Matrix m1({{12,2, 4, 5}, {3, 2, 1, 1}, {1, 1, 5, 6}});
        QRDecomposer qr4(m1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    if (brojac==3) std::cout << "Test QRDecomposer: OK" << std::endl;
    else std::cout << "Test QRDecomposer: Greška" << std::endl;
}

void QRSolveTest() {
    int brojac(0);
    Matrix A({{2, 4}, {3, -3}});
    QRDecomposer qr(A);
    Vector b({-2, 5});
    Vector x(2);
    qr.Solve(b, x);
    if (x.EqualTo({0.777778, -0.888889}, 0.0001)) brojac++;
    if (qr.Solve(b).EqualTo({0.777778, -0.888889}, 0.0001)) brojac++;
    
    //Incompatible formats:
    Vector b2({1, 2, 3});
    try {
        qr.Solve(b2, x);
    }
    catch(std::domain_error) {brojac++;}
    catch(...) {}
    
    //LUSolve sa matricama:
    Matrix m1({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Matrix b1({{4, 4, 4}, {8, 7, 1}, {2, 2, 3}});
    Matrix x1(3,3);
    QRDecomposer qr1(m1);
    qr1.Solve(b1, x1);
    
    Matrix xrez({{5.08108, 3.59459, -5.43243}, {0.27027, 0.648649, 2.89189}, {-0.540541, -0.297297, 1.21622}});
    if (x1.EqualTo(xrez, 0.0001)) brojac++;
    if (qr1.Solve(b1).EqualTo(xrez, 0.0001)) brojac++;
    
    if (brojac==5) std::cout << "Test Solve QR: OK" << std::endl;
    else std::cout << "Test Solve QR: Greška" << std::endl;
}

void MulQWithTest() {
    Matrix m({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Vector x{1, -3, 2};
    QRDecomposer qr(m);
    Matrix Q=qr.MulQWith(x);
    if (Q.EqualTo(qr.GetQ()*x)) std::cout << "Test Mul Q With: OK" << std::endl;
    else std::cout << "Test Mul Q With: Greška" << std::endl;
}

void MulQTWithTest() {
    Matrix m({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Vector x{1, -3, 2};
    QRDecomposer qr(m);
    Matrix Q=qr.MulQTWith(x);
    if (Q.EqualTo(Transpose(qr.GetQ())* x)) std::cout << "Test Mul QT With: OK" << std::endl;
    else std::cout << "Test Mul QT With: Greška" << std::endl;
}



int main ()
{
    VectorChopTest();
    VectorEqualToTest();
    
    
    MatrixChopTest();
    MatrixEqualToTesst();
    LeftDivTest();
    RightDivTest();
    DetTest();
    InvertTest();
    RREFTest();
    
    
    LUDecomposerTest();
    SolveLUTest();
    
    QRDecomposerTest();
    QRSolveTest();
    MulQTWithTest();
    MulQTWithTest();
	return 0;
}

