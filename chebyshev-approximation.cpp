//NA 2017/2018: Zadaća 4, Zadatak 1
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <utility>
#include <limits>


const double pi=4*atan(1);

template <typename FunType>
std::pair<double, bool> RombergIntegration(FunType f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 50) {
    if (eps<0 || nmax<0 || nmin<0 || nmax<nmin) throw std::domain_error("Bad parameter");
    
    int N(2);
    double h((b-a)/N), s((f(a)+f(b))/2), p(0);
    double iold(s);
    std::vector<double> in; 
    
    for (int i=1; N<=nmax; i++) {
        for (int j=1; j<=N/2; j++) s=s+f(a+(2*j-1)*h);
        in.push_back(h*s);
        p=4;
        for (int k=in.size()-2; k>=0; k--) {
            in[k]=(p*in[k+1]-in[k])/(p-1);
            p=p*4;
        }
        if (std::fabs(in[0]-iold)<=eps) return std::make_pair(in[0], true);
        iold=in[0];
        h=h/2;
        N=N*2;
    }
    return std::make_pair(iold, false);
}

template <typename FunType>
std::pair<double, bool> TanhSinhIntegration(FunType f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 20, double range = 3.5) {
    if (eps<0 || nmax<0 || nmin<0 || nmax<nmin || range<0) throw std::domain_error("Bad parameter");
    int N(2);
    double h(2*range/N), p((b+a)/2), q((b-a)/2), s(0), t(0), u(0), v(0);
    double iold(s), in(0);
    
    while (N<nmax) {
        for (int i=1; i<=N/2; i++) {
            t=-range+(2*i-1)*h;
            u=(pi*std::sinh(t))/2;
            v=f(p+q*std::tanh(u));
            if (std::isfinite(v)) s=s+(q*pi*std::cosh(t)*v)/(2*std::cosh(u)*std::cosh(u));
        }
        in=h*s;
        if (N>=nmin && std::fabs(in-iold)<=eps) return std::make_pair(in, true);
        iold=in;
        N=2*N;
        h=h/2;
    }
    return std::make_pair(in, false);
}


template <typename FunType>
std::pair<double,bool> AdaptiveAux(FunType f, double a, double b, double eps, double f1, double f2, double f3, int maxdepth) {
    if(!std::isfinite(f1)) f1=0;
    if(!std::isfinite(f2)) f2=0;
    if(!std::isfinite(f3)) f3=0;
    
    double c=(a+b)/2;
    double i1=(b-a)*(f1+4*f3+f2)/6;
    double f4=f((a+c)/2), f5=f((c+b)/2);
    
    if(!std::isfinite(f4)) f4=0;
    if(!std::isfinite(f5)) f5=0;
    
    double i2=((b-a)*(f1+4*f4+2*f3+4*f5+f2))/12;
    
    if(std::abs(i1-i2)<=eps) return std::make_pair(i2, true);
    if(maxdepth<=0) return std::make_pair(i2, false);
    
    auto aa1=AdaptiveAux(f, a, c, eps, f1, f3, f4, maxdepth-1);
    auto aa2=AdaptiveAux(f, c, b, eps, f3, f2, f5, maxdepth-1);
    
    return std::make_pair(aa1.first+aa2.first, aa1.second && aa2.second);
}

template <typename FunType>
std::pair<double, bool> AdaptiveIntegration(FunType f, double a, double b, double eps = 1e-10, int maxdepth = 30, int nmin = 1) {
    if(eps<0 || maxdepth<0 || nmin<0) throw std::domain_error("Bad parameter");
    
    std::pair<double,bool> s(std::make_pair(0, true));
    double h=(b-a)/nmin;
    for(int i=1; i<=nmin; i++) {
        auto aa=AdaptiveAux(f, a, a+h, eps, f(a), f(a+h), f(a+h/2), maxdepth);
        s=std::make_pair( s.first+aa.first, s.second && aa.second );
        a=a+h;
    }
    return s;
}


class ChebyshevApproximation {
    std::vector<double> c;
    double xmin, xmax;
    int m;
    int n;
    
    ChebyshevApproximation(std::vector<double> c1, double xmin, double xmax, int n) : c(c1), xmin(xmin), xmax(xmax), m(c1.size()-1), n(n) {} 
    
    public:
    template <typename FunType>
    ChebyshevApproximation(FunType f, double xmin, double xmax, int n) : xmin(xmin), xmax(xmax), m(n), n(n) {
        if (xmin>=xmax || n<1) throw std::domain_error("Bad parameters");

        std::vector<double> w(n+2);
        std::vector<double> v(n+1);
        
        for (int i=0; i<=n+1; i++) w[i]=std::cos((pi*i)/(2*n+2));
        for (int i=0; i<=n/2; i++) v[i]=f((xmin+xmax+(xmax-xmin)*w[2*i+1])/2);
        for (int i=n/2+1; i<=n; i++) v[i]=f((xmin+xmax-(xmax-xmin)*w[2*n+1-2*i])/2);
        
        
        c.resize(n+1);
        double s(0), p(0);
        for (int k=0; k<=n; k++) {
            s=0;
            for (int i=0; i<=n; i++) {
                p=(k*(2*i+1))%(4*n+4);
                if (p>2*n+2) p=4*n+4-p;
                if (p>n+1) s=s-v[i]*w[2*n+2-p];
                else s=s+v[i]*w[p];
            }
            c[k]=2*s/(n+1);
        }
        
    }
    
    void set_m(int m) {
        if (m<=1 || m>c.size()-1) throw std::domain_error("Bad order");
        this->m=m;
    }
    void trunc(double eps) {
        if (eps<0) throw std::domain_error("Bad tolerance");
        int i;
        for (i=m; i>=0; i--) if (std::fabs(c[i])>=eps) break;
        if (i==-1) throw std::domain_error("Bad tolerance");
        m=i;
    }
    double operator()(double x) const {
        if (x<xmin || x>xmax) throw std::domain_error("Bad argument");
        
        double t, p(0), q(c[m]), r(0);
        t=(2*x-xmin-xmax)/(xmax-xmin);
        for (int k=m-1; k>=1; k--) {
            r=2*t*q-p+c[k];
            p=q;
            q=r;
        }
        return (t*q-p+c[0]/2);
        
    }
    
    double derivative(double x) const {
        if (x<xmin || x>xmax) throw std::domain_error("Bad argument");
        double t((2*x-xmin-xmax)/(xmax-xmin)), p(1), q(4*t), s(c[1]+4*c[2]*t), r(0);
        for(int k=3; k<=m; k++) {
            r=k*((2*t*q)/(k-1)-p/(k-2));
            s=s+c[k]*r;
            p=q;
            q=r;
        }
        return 2*s/(xmax-xmin);
    }

    ChebyshevApproximation derivative() const {
        std::vector<double> c1(c.size());
        double mi(4./(xmax-xmin));
        c1[m-1]=mi*m*c[m];
        c1[m-2]=mi*(m-1)*c[m-1];
        for (int k=m-3; k>=0; k--) c1[k]=c1[k+2]+mi*(k+1)*c[k+1];
        return ChebyshevApproximation(c1, xmin, xmax, n);
    }
    
    ChebyshevApproximation antiderivative() const {
        std::vector<double> c1(m+2);
        double mi((xmax-xmin)/4);
        c1[0]=0;
        for (int i=1; i<=m-1; i++) c1[i]=mi/i*(c[i-1]-c[i+1]);
        c1[m]=mi/m*c[m-1];
        c1[m+1]=mi/(m+1)*c[m];
        return ChebyshevApproximation(c1, xmin, xmax, m+2);
    }
    
    double integrate(double a, double b) const {
        if (a<xmin || a>xmax || b<xmin || b>xmax) throw std::domain_error("Bad interval");
        
        ChebyshevApproximation F(antiderivative());
        return (F(b)-F(a));
    }
    
    double integrate() const {
        double s(0);
        for (int i=1; i<=(m-1)/2+1; i++) s=s+2*c[2*i]/(1-4*i*i);
        s=s*(xmax-xmin)/2;
        s=s+(c[0]*(xmax-xmin))/2;
        return s;
    }
};

bool doubleequal(double x, double y) {
    double eps=10*std::numeric_limits<double>::epsilon()*(std::fabs(x)+std::fabs(y));
    return (std::fabs(x-y)<=eps && !((x<0 && y>0) || (x>0 && y<0)));
} 


void testChebyshevApproximation() {
    int brojac(0);
    
    auto f=[](double x){return x*x*x;};
    
    //Bad parameters
    try {
        ChebyshevApproximation(f, 20, 10, 10);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    auto cheby=ChebyshevApproximation(f, 0, 10, 20);
    cheby.trunc(0.001);
    
    //Bad argument
    try {
        cheby(11);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    if (doubleequal(27, cheby(3))) brojac++;
    if (doubleequal(cheby.derivative(2), 12)) brojac++;
    if (doubleequal(cheby.integrate(2, 2), 0)) brojac++;
    if (doubleequal(cheby.derivative(1.7), cheby.derivative()(1.7))) brojac++;
    
    if (brojac==6) std::cout << "Test ChebyshevApproximation: OK" << std::endl;
    else std::cout << "Test ChebyshevApproximation: Greska" << std::endl;
}


void testRombergIntegration() {
    int brojac(0);
    
    auto f=[](double x){ return std::sin(x); };
    auto in=[](double x){ return -std::cos(x); };
    
    //Bad parameter
    try {
        RombergIntegration(f, 0, pi, -0.1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    
    auto rez = RombergIntegration(f, 1.5708, 0, 1e-8, 10000);
    
    if (std::fabs(rez.first-in(0))<0.0001) brojac++;
    if (rez.second==true) brojac++;
    
    if (brojac==3) std::cout << "Test RombergIntegration: OK" << std::endl;
    else std::cout << "Test RombergIntegration: Greska" << std::endl;
}

void testTanhSinhIntegration() { 
    int brojac(0);
    
    auto f=[](double x){ return std::sin(x)*std::sin(x)*std::cos(x); };
    auto in=[](double x){ return std::sin(x)*std::sin(x)*std::sin(x)/3; };
    
    //Bad parameter
    try {
        RombergIntegration(f, 0, pi, -0.1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    auto rez = RombergIntegration(f, 0, 0.5, 1e-8, 10000);
    
    if (std::fabs(rez.first-in(0.5))<0.0001) brojac++;
    if (rez.second==true) brojac++;
    
    if (brojac==3) std::cout << "Test TanhSinhIntegration: OK" << std::endl;
    else std::cout << "Test TanhSinhIntegration: Greska" << std::endl;
}


void testAdaptiveIntegration() {
    int brojac(0);
    
    auto f=[](double x){ return 1./std::sqrt(x); };
    auto in=[](double x){ return 2*std::sqrt(x); };
    
    //Bad parameter
    try {
        AdaptiveIntegration(f, 0, 1, -0.1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    auto rez = AdaptiveIntegration(f, 0, 1, 1e-8, 10000);
    
    if (std::fabs(rez.first-in(1))<0.0001) brojac++;
    if (rez.second==true) brojac++;
    
    if (brojac==3) std::cout << "Test AdaptiveIntegration: OK" << std::endl;
    else std::cout << "Test AdaptiveIntegration: Greska" << std::endl;
}

int main ()
{
    testChebyshevApproximation();
    testRombergIntegration();
    testTanhSinhIntegration();
    testAdaptiveIntegration();
    return 0;
}