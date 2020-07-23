//NA 2018/2019: Zadaća 3, Zadatak 1
#include <iostream>
#include <cmath>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <vector>
#include <limits>

typedef std::pair<double, double> par;

bool doubleequal(double x, double y) {
    double eps=10*std::numeric_limits<double>::epsilon()*(std::fabs(x)+std::fabs(y));
    return (std::fabs(x-y)<=eps && !((x<0 && y>0) || (x>0 && y<0)));
} 


class AbstractInterpolator {
    protected:
    std::vector<par> v;
    mutable int cache=1;
    
    int Locate(double x) const {
        
        if (x<=v.at(0).first || doubleequal(x, v.at(0).first)) return 0; // x <= y
        if (x>v.at(v.size()-1).first && !doubleequal(x, v.at(v.size()-1).first)) return v.size(); // x > y
        
        
        if (x>v.at(cache-1).first && x<=v.at(cache).first) return cache;                             //kao i prosli
        if (cache>1 && v.size()>2) {
            if (x>v.at(cache-2).first && x<=v.at(cache-1).first) { cache--; return cache; }          //prvi prethodni
        }
        if (cache<v.size()-1 && v.size()>2) {
            if (x>v.at(cache).first && x<=v.at(cache+1).first) { cache++; return cache; }            //prvi sljedeci
        }
        
        
        auto low=std::lower_bound(v.begin(), v.end(), x, [](par p, double target){ return target>p.first && !doubleequal(target, p.first); });
        cache=low-v.begin();
        return cache;
    }
    
    public:
    AbstractInterpolator(const std::vector<par> &data) {
        if (data.size()==0 ||  data.size()==1) throw std::domain_error("Invalid data set");
        v=data;
        std::sort(v.begin(), v.end(), [](par p1, par p2){ return p1.first<p2.first; } );
        if (std::adjacent_find(v.begin(), v.end(), [](par p1, par p2){ return doubleequal(p1.first, p2.first); })!=v.end()) throw std::domain_error("Invalid data set");
    }
    
    virtual double operator()(double x) const = 0;
};


class LinearInterpolator : public AbstractInterpolator {
    
    public:
    LinearInterpolator(const std::vector<par> &data) : AbstractInterpolator(data) {}
    double operator()(double x) const override {
        int index=Locate(x);
        
        if (index==0) index=1;
        else if (index==v.size()) index=v.size()-1;
        
        double y=(v.at(index).first-x)*v.at(index-1).second/(v.at(index).first-v.at(index-1).first)+((x-v.at(index-1).first)*v.at(index).second)/(v.at(index).first-v.at(index-1).first);
        return y;
    }
};

void testLinearInterpolator() {
    int brojac(0);
    LinearInterpolator LI({ {1, 2}, {7, 4}, {3, 3}, {9, 12}, {8, 10} });
    
    if (doubleequal(LI(3.5), 3.125)) brojac++;
    
    //Rubne tacke:
    if (doubleequal(LI(7),4)) brojac++;
    if (doubleequal(LI(9), 12)) brojac++;
    
    //Ekstrapolacija
    if (doubleequal(LI(-5), -1)) brojac++;
    if (doubleequal(LI(18), 30)) brojac++;
    
    //Iste tacke:
    try {
        LinearInterpolator LI2({ {1, 2}, {2, 7}, {8, 3}, {2, 6}, {12, 11} });
    }
    catch(std::domain_error) {
        brojac++;
    }
    catch(...) {}
    
    if (brojac==6) std::cout << "Test LinearInterpolator: OK" << std::endl;
    else std::cout << "Test LinearInterpolator: Greska" << std::endl;
}


class PolynomialInterpolator : public AbstractInterpolator {
    std::vector<double> y;
    
    public:
    PolynomialInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data) {
        
        y.push_back(v.at(v.size()-1).second);
        
        for (int j=1; j<=v.size()-1; j++) {
            for (int i=v.size(); i>=j+1; i--) {
                v.at(i-1).second=(v.at(i-1).second-v.at(i-2).second)/(v.at(i-1).first-v.at(i-j-1).first);
            }
            y.push_back(v.at(v.size()-1).second);
        }
    }
    
    double operator()(double x) const override {
        
        double f(v.at(v.size()-1).second);
        for (int i=v.size()-1; i>=1; i--) {
            f=f*(x-v.at(i-1).first)+v.at(i-1).second;
        }
        return f;
        
    }
    
    void AddPoint(const std::pair<double, double> &p) {
        if (std::any_of(v.begin(), v.end(), [p](par p1) { return doubleequal(p1.first, p.first); })) throw std::domain_error("Invalid point");
        
        v.push_back(p);
        
        int brEl(v.size());
        y.resize(brEl);
        for (int i=1; i<=brEl-1; i++) {
            double Yi=y.at(i-1);
            y.at(i-1)=v.at(brEl-1).second;
            v.at(brEl-1).second=(v.at(brEl-1).second-Yi)/(v.at(brEl-1).first-v.at(brEl-i-1).first);
        }
        y.at(brEl-1)=v.at(brEl-1).second; 
        
    }
    
    std::vector<double> GetCoefficients() const {
        
        int n(v.size());
        
        std::vector<double> master(n+1),polinom(n+1);
        master.at(0)=1;
        for(int i=1; i<=n; i++) {
            for(int j=0; j<=i; j++)
                polinom.at(j)+=v.at(i-1).second*master.at(j);
            master.at(i)=master.at(i-1);
            for(int j=i-1; j>0; j--) master.at(j)=master.at(j-1)-v.at(i-1).first*master.at(j);
            master.at(0)*=(-1)*v.at(i-1).first;    
        }
        polinom.resize(n);
        return polinom;
    }
};


void testPolynomialInterpolator() {
    int brojac(0);
    auto f([](double x){ return 5*x*x+6*x+12; } );
    
    PolynomialInterpolator PI({ {1, f(1)}, {3, f(3)}, {8, f(8)} });
    
    //Operator ():
    
    if (doubleequal(PI(2), 44)) brojac++;
    if (doubleequal(PI(-1), 11)) brojac++;
    
    
    //GetCoefficients:
    auto v(PI.GetCoefficients());
    std::vector<double> koef={12, 6, 5};
    if (v==koef) brojac++;
    
    //AddPoint
    PI.AddPoint({11, f(11)});
    koef={12, 6, 5, 0};
    auto v1(PI.GetCoefficients());
    if (v1==koef) brojac++;
    
    //AddPoint izuzetak:
    try {
        PI.AddPoint({3, f(3)});
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    if (brojac==5) std::cout << "Test PolynomialInterpolator: OK" << std::endl;
    else std::cout << "Test PolynomialInterpolator: Greska" << std::endl;
}


class PiecewisePolynomialInterpolator : public AbstractInterpolator {
    int k;
    
    public:
    PiecewisePolynomialInterpolator (const std::vector<par> &data, int order) : AbstractInterpolator(data) {
        if (order<1 || order>=data.size()) throw std::domain_error("Invalid order");
        k=order;
    }
    
    double operator()(double x) const override {
        int index=Locate(x);
        
        if (index==0) index=1;
        else if (index==v.size()) index=v.size()-1;
        
        int beg, end;
        
        if (k%2==0) {
            beg=index-k/2;
            end=index+k/2;
        } 
        else {
            beg=index-(k-1)/2;
            end=index+(k+1)/2;
        }
        
        if (beg<1) { beg=1; end=k+1; }
        else if (end>v.size()) { end=v.size(); beg=end-k; }
        
        double sum=0, prod=1;
        for (int i=beg; i<=end; i++) {
            prod=1;
            for (int j=beg; j<=end; j++) {
                if (j!=i) prod*=(x-v.at(j-1).first)/(v.at(i-1).first-v.at(j-1).first);
            }
            sum+=v.at(i-1).second*prod;
        }
        return sum;
    }
};

void testPiecewisePolynomialInterpolator() {
    
    auto f([](double x){ return 1./(1+x*x); });
    
    std::vector<par> data({ {0,f(0)}, {1,f(1)}, {2,f(2)}, {3,f(3)}, {4,f(4)}, {5, f(5)}, {6, f(6)} });
    int brojac(0);
    
    //Red < 1
    try {
        PiecewisePolynomialInterpolator PPI(data, 0);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    //Red > broj tacaka
    try {
        PiecewisePolynomialInterpolator PPI(data, 8);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    PiecewisePolynomialInterpolator PPI(data, 3);
    
    //Operator()
    if (std::fabs(PPI(5.3)-f(5.3))<0.01) brojac++;
    if (doubleequal(PPI(3), f(3))) brojac++;
    if (std::fabs(PPI(6.1)-f(6.1))<0.01) brojac++;
    
    if (brojac==5) std::cout << "Test PiecewisePolynomialInterpolator: OK" << std::endl;
    else std::cout << "Test PiecewisePolynomialInterpolator: Greska" << std::endl;
}


class SplineInterpolator : public AbstractInterpolator {
    std::vector<double> r;
    
    public:
    SplineInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data) {
        int n(v.size());
        
        r.resize(n, 0);
        
        std::vector<double> a(n-2);
        
        for (int i=1; i<n-1; i++) {
            a.at(i-1)=2*(v.at(i+1).first-v.at(i-1).first);
            r.at(i)=3*((v.at(i+1).second-v.at(i).second)/(v.at(i+1).first-v.at(i).first)-(v.at(i).second-v.at(i-1).second)/(v.at(i).first-v.at(i-1).first));
        }
        
        for (int i=1; i<n-2; i++) {
            double mi((v.at(i+1).first-v.at(i).first)/a.at(i-1));
            a.at(i)=a.at(i)-mi*(v.at(i+1).first-v.at(i).first);
            r.at(i+1)=r.at(i+1)-mi*r.at(i);
        }
        
        r.at(n-2)/=a.at(n-3);
        
        for (int i=n-3; i>=1; i--) {
            r.at(i)=(r.at(i)-(v.at(i+1).first-v.at(i).first)*r.at(i+1))/a.at(i-1);
        }
    }
    
    double operator()(double x) const override {
        int index=Locate(x);
        
        if (index==0) index=0;
        else if (index==v.size()) index=v.size()-2;
        else index--;
        
        double s((r.at(index+1)-r.at(index))/(3*(v.at(index+1).first-v.at(index).first)));
        double q( (v.at(index+1).second-v.at(index).second)/(v.at(index+1).first-v.at(index).first) - (v.at(index+1).first-v.at(index).first)*(r.at(index+1)+2*r.at(index))/3. );
        //pi=yi
        
        double X (x - v.at(index).first);
        
        return v.at(index).second + X * ( q + X * (r.at(index) + s*X));
    }

};

void testSplineInterpolator() {
    
    double pi(4*atan(1));
    int brojac(0);
    auto f([](double x){ return std::sin(x); });
    
    std::vector<par> data({ {0, f(0)}, {pi/2, f(pi/2)}, {pi, f(pi)}, {3*pi/2, f(3*pi/2)}, {2*pi, f(2*pi)}, {5*pi/2, f(5*pi/2)}, {3*pi, f(3*pi)} });
    
    SplineInterpolator SI(data);
    
    if (std::fabs(f(3)-SI(3))<0.01) brojac++;
    if (std::fabs(f(1.7)-SI(1.7))<0.01) brojac++;
    if (std::fabs(f(-1.7)-SI(-1.7))<0.01) brojac++;
    if (std::fabs(f(4.99)-SI(4.99))<0.01) brojac++;
    
    if (brojac==4) std::cout << "Test SplineInterpolator: OK" << std::endl;
    else std::cout << "Test SplineInterpolator: Greska" << std::endl;
}


class BarycentricInterpolator : public AbstractInterpolator {
    int d;
    std::vector<double> w;
    
    public:
    BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order) : AbstractInterpolator(data) {
        if (order<0 || order>data.size()) throw std::domain_error("Invalid order");
        d=order;
        int n(v.size());
        w.resize(n);
        double p(1);
        
        for (int i=0; i<n; i++) {
            w.at(i)=0;
            for (int k=std::max(0, i-d); k<=std::min(i, n-d-1); k++) {
                p=1;
                for (int j=k; j<k+d; j++) {
                    if (i!=j) p/=(v.at(i).first-v.at(j).first);
                }
                if (k%2==1) p*=-1;
            }
            w.at(i)+=p;
        }
    }
    
    double operator()(double x) const override {
        double p(0), q(0), u(0);
        
        for (int i=0; i<v.size(); i++) {
            if (doubleequal(x, v.at(i).first)) return v.at(i).second;
            u=w.at(i)/(x-v.at(i).first);
            p+=u*v.at(i).second;
            q+=u;
        }
        return p/q;
    }
    
    std::vector<double> GetWeights() const { return w; }
};

void testBarycentricInterpolator() {
    
    auto f([](double x){ return std::sin(x)*std::sin(x); });
    int brojac(0);
    double pi(4*atan(1));
    
    std::vector<par> data({ {0, f(0)}, {pi/4, f(pi/4)}, {pi/2, f(pi/2)}, {3*pi/4, f(3*pi/4)}, {pi, f(pi)}, {5*pi/4, f(5*pi/4)}, {3*pi/2, f(3*pi/2)}, {7*pi/4, f(7*pi/4)}, {2*pi, f(2*pi)} });
    
    //Red < 0
    try {
        BarycentricInterpolator BI(data, -1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    //Red > broj tacaka
    try {
        BarycentricInterpolator BI(data, 10);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    BarycentricInterpolator BI(data, 0);
    
    //GetWeights
    std::vector<double> w({ 1, -1, 1, -1, 1, -1, 1, -1, 1 });
    if (BI.GetWeights()==w) brojac++;
    
    //Operator ()
    if (std::fabs(f(0.5)-BI(0.5))<0.01) brojac++;
    if (std::fabs(f(2.6)-BI(2.6))<0.01) brojac++;
    if (std::fabs(f(3.4)-BI(3.4))<0.01) brojac++;
    if (std::fabs(f(3.88)-BI(3.88))<0.01) brojac++;
    if (std::fabs(f(5.7)-BI(5.7))<0.01) brojac++;
    
    if (brojac==8) std::cout << "Test BarycentricInterpolator: OK" << std::endl;
    else std::cout << "Test BarycentricInterpolator: Greska" << std::endl;
}


template <typename FunType>
std::pair<double, bool> Limit (FunType f, double x0, double h = 0, double eps = 1e-8, double nmax = 20) {
    if (eps<=0 || nmax<3 || nmax>30) throw std::domain_error("Invalid parameters");
    
    bool funkcija(false);
    
    auto g=[f](double x){ return f(1./x); };
    
    if (!std::isinf(x0) && doubleequal(h, 0)) h=0.001*std::max(1., std::fabs(x0));
    
    if (std::isinf(x0)) {
        if (x0<0) x0=-std::numeric_limits<double>::min();
        else x0=0;
        funkcija=true;
    }
    
    double p;
    std::vector<double> y(nmax);
    double yold(std::numeric_limits<double>::infinity());
    
    for (int i=0; i<nmax; i++) {
        if (funkcija) y.at(i)=g(x0+h);
        else y.at(i)=f(x0+h);
        p=2;
        for (int k=i-1; k>=0; k--) {
            y.at(k)=(p*y.at(k+1)-y.at(k))/(p-1);
            p=2*p;
        }
        if (std::fabs(y.at(0)-yold)<eps) return std::make_pair(y.at(0), true);
        yold=y.at(0);
        h=h/2;
    }
    return std::make_pair(y.at(0), false);
}


void testLimit() {
    int brojac(0);
    
    //Neispravni parametri:
    try {
        Limit([](double x){ return std::sin(x); } , 0, 0, -1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    try {
        Limit([](double x){ return std::sin(x); } , 0, 0, 0.0001, 2);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    try {
        Limit([](double x){ return std::sin(x); } , 0, 0, 0.0001, 31);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    double inf = std::numeric_limits<double>::infinity();
    
    //Neki tablicni limesi:
    auto limes =  Limit([](double x) { return std::sin(x); }, 0);
    if (std::fabs(limes.first)<0.00001 && limes.second==true) brojac++;
    
    limes =  Limit([](double x) { return std::cos(x); }, 0);
    if (std::fabs(limes.first-1)<0.00001 && limes.second==true) brojac++;
    
    limes =  Limit([](double x) { return std::sin(x)/x; }, 0);
    if (std::fabs(limes.first-1)<0.00001 && limes.second==true) brojac++;
    
    limes =  Limit([](double x) { return std::tan(x)/x; }, 0);
    if (std::fabs(limes.first-1)<0.00001 && limes.second==true) brojac++;
    
    limes =  Limit([](double x) { return x/(x*x); }, inf, -1);
    if (std::fabs(limes.first)<0.00001 && limes.second==true) brojac++;
    
    limes =  Limit([](double x) { return std::pow(1+1./x, x); }, inf, 1);
    if (std::fabs(limes.first-std::exp(1))<0.00001 && limes.second==true) brojac++;
    
    limes =  Limit([](double x) { return x*x*x/x; }, -inf);
    if (limes.second==false) brojac++;
    
    limes =  Limit([](double x) { return (1-std::cos(x))/(x*x); }, 0);
    if (std::fabs(limes.first-0.5)<0.00001 && limes.second==true) brojac++;
    
    if (brojac==11) std::cout << "Test Limit: OK" << std::endl;
    else std::cout << "Test Limit: Greska" << std::endl;
}



int main ()
{
    testLinearInterpolator();
    testPolynomialInterpolator();
    testPiecewisePolynomialInterpolator();
	testSplineInterpolator();
	testBarycentricInterpolator();
	testLimit();
	return 0;
}
