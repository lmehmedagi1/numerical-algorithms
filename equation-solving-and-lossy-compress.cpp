//NA 2017/2018: Zadaća 5, Zadatak 1
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>
#include <complex>
#include <utility>
#include <algorithm>

#include <iomanip>

using namespace std;


const double pi=4*std::atan(1);

enum RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic};

bool doubleequal(double x, double y) {
    double eps=10*std::numeric_limits<double>::epsilon()*(std::fabs(x)+std::fabs(y));
    return (std::fabs(x-y)<=eps && !((x<0 && y>0) || (x>0 && y<0)));
} 

int sgn(double val) {
    if (doubleequal(val, 0)) return 0;
    return (0 < val) - (val < 0);
}

double randBr (double min, double max) {
    double x = (double)rand() / RAND_MAX;
    return min + x * (max - min);
}
std::complex<double> RandomComplex(double min1, double max1, double min2, double max2) {
    return { randBr(min2, max2), randBr(min2, max2) };
}

std::pair<std::complex<double>, bool> Laguerre(std::vector<std::complex<double>> p, int n, std::complex<double> x, double eps, int maxiter) {
    std::complex<double> deltax=std::numeric_limits<double>::infinity();
    int k(1);
    std::complex<double> f, d, s, r;
    
    while (std::abs(deltax)>eps && k<maxiter) {
        f=p.at(n);
        d=0;
        s=0;
        for (int i=n-1; i>=0; i--) {
            s=s*x+2.*d;
            d=d*x+f;
            f=f*x+p.at(i);
        }
        if (std::abs(f)<=eps) return {x, true};
        r=std::sqrt(double(n-1)*(double(n-1)*d*d-double(n)*f*s));
        
        if (std::abs(d+r)>std::abs(d-r)) deltax=double(n)*f/(d+r);
        else deltax=double(n)*f/(d-r);
        
        x=x-deltax;
        k++;
    }
    if (std::abs(deltax)<=eps) return {x, true};
    return {x, false};
}
std::pair<std::complex<double>, bool> Laguerre(std::vector<double> p, int n, std::complex<double> x, double eps, int maxiter) {
    std::vector<std::complex<double>> p1(p.size());
	for (int i=0; i<p1.size(); i++) p1[i] = std::complex<double>(p[i]);
	return Laguerre(p1, n, x, eps, maxiter);
}

bool isti(std::vector<std::complex<double>> v1, std::vector<std::complex<double>> v2) {
    if (v1.size()!=v2.size()) return false;
    for (int i=0; i<v1.size(); i++) {
        if (std::abs(v1[i]-v2[i])>0.00001) return false;
    }
    return true;
}

template <typename FunTip>
double RK4Step(FunTip f, double x, double y, double h) {
    double K1=f(x, y);
	double K2=f(x+h/2, y+h*K1/2);
	double K3=f(x+h/2, y+h*K2/2);
	double K4=f(x+h, y+h*K3);
	return y+h*(K1+2*K2+2*K3+K4)/6;
}

template <typename FunTip>
bool BracketRoot(FunTip f, double x0, double &a, double &b, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4) {
    if (hinit <= 0 || hmax <= 0 || lambda <= 0) throw std::domain_error("Invalid parameters");
    
    bool pronadjen(true);
    
    double a1(x0),b1;
    
    double h = hinit;
    double f1(f(a1));
    double f2(0);
    
    while (std::fabs(h) < hmax) {
        b1 = a1 + h;
        f2 = f(b1);
        while (!std::isfinite(f2)) {
            h = h/(2*(1+lambda));
            if (std::fabs(h)<=std::fabs(a1)*std::numeric_limits<double>::epsilon()) {
                pronadjen = false;
                break;
            }
            b1=a1+h;
            f2=f(b1);
        }
        if (!pronadjen) break;
        if (f1*f2<=0) {
            if (a1>b1) std::swap(a1, b1);
            a=a1;
            b=b1;
            return true;
        }
        h=lambda*h;
        a1=b1;
        f1=f2;
    }
    
    a1 = x0;
    h = -hinit;
    f1 = f(a1);
    
    while (std::fabs(h) < hmax) {
        b1 = a1 + h;
        f2 = f(b1);
        while (!std::isfinite(f2)) {
            h = h/(2*(1+lambda));
            if (std::fabs(h)<=std::fabs(a1)*std::numeric_limits<double>::epsilon()) return false; 
            b1=a1+h;
            f2=f(b1);
        }
        if (f1*f2<=0) {
            if (a1>b1) std::swap(a1, b1);
            a=a1;
            b=b1;
            return true;
        }
        h=lambda*h;
        a1=b1;
        f1=f2;
    }
    
    return false;
}

template <typename FunTip>
double RegulaFalsiSolve(FunTip f, double a, double b, RegulaFalsiMode mode = Slavic, double eps = 1e-10, int maxiter = 100) {
    if (f(a)*f(b)>0) throw std::range_error("Root must be bracketed");
    if (maxiter<=0 || eps<=0) throw std::domain_error("Invalid parameters");
    
    double f1, f2, f3(0), c(a), cold(b);
    
    auto phi = [](double x){return x/(1+std::fabs(x));};
    
    if (mode==Unmodified || mode==Illinois) { f1=f(a); f2=f(b); }
    if (mode==Slavic || mode==IllinoisSlavic) { f1=phi(f(a)); f2=phi(f(b)); }
    
    
    for (int i=0; i<maxiter; i++) {
        if (std::fabs(c-cold) > eps) {
            cold=c;
            c=(a*f2 - b*f1)/(f2-f1);
            
            if (mode==Unmodified || mode==Illinois) f3=f(c);
            else if (mode==Slavic || mode==IllinoisSlavic) f3=phi(f(c));
            
            if (doubleequal(f3, 0)) return c;
            if (f1*f3<0) {
                b=a;
                f2=f1;
            }
            else if (mode==Illinois || mode==IllinoisSlavic) {
                f2=f2/2;
            }
            a=c;
            f1=f3;
        }
        else return c;
    }
    
    throw std::logic_error("Given accuracy has not achieved");
}

template <typename FunTip>
double RiddersSolve(FunTip f, double a, double b, double eps = 1e-10, int maxiter = 100) {
    if (f(a)*f(b)>0) throw std::range_error("Root must be bracketed");
    if (maxiter<=0 || eps<=0) throw std::domain_error("Invalid parameters");
    
    
    double f1(f(a)), f2(f(b)), f3(0), f4(0), c(0), d(0);
    
    for (int i=0; i<maxiter; i++) {
        if (std::fabs(b-a) > eps) {
            c=(a+b)/2;
            f3=f(c);
            if (doubleequal(f3, 0)) return c;
            d=c+(f3*(c-a)*sgn(f1-f2))/std::sqrt(f3*f3-f1*f2);
            f4=f(d);
            if (doubleequal(f4, 0)) return d;
            if (f3*f4<=0) {
                a=c;
                b=d;
                f1=f3;
                f2=f4;
            }
            else if (f1*f4<=0) {
                b=d;
                f2=f4;
            }
            else {
                a=d;
                f1=f4;
            }
        }
        else return (a+b)/2;
    }
    
    throw std::logic_error("Given accuracy has not achieved");
    
} 

template <typename FunTip1, typename FunTip2>
double NewtonRaphsonSolve(FunTip1 f, FunTip2 fprim, double x0, double eps = 1e-10, double damping = 0, int maxiter = 100) {
    
    if (maxiter<=0 || eps<=0 || damping<0 || damping>=1) throw std::domain_error("Invalid parameters");
    
    
    double deltax = std::numeric_limits<double>::infinity();
    double x = x0;
    double v(f(x)), d(fprim(x)), w(0);
    
    for (int i=0; i<maxiter; i++) {
        if (std::fabs(deltax)>eps) {
            
            if (std::fabs(v)<=eps) return x;
            deltax = v/d;
            
            w = v;
            v=f(x-deltax);
            d=fprim(x-deltax);
            
            while (std::fabs(v) > std::fabs(w) || !std::isfinite(v) || doubleequal(d, 0)) {
                if (doubleequal(damping, 0)) throw std::logic_error("Convergence has not achieved");
                deltax=damping*deltax;
                v=f(x-deltax);
                d=fprim(x-deltax);
            }
            x=x-deltax;
        }
        else return x;
    }
    throw std::logic_error("Convergence has not achieved");
}

std::vector<std::complex<double>>PolyRoots(std::vector<std::complex<double>> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10) {
    if (eps<=0 || maxiters<=0 || maxtrials<=0) throw std::domain_error("Invalid parameters");
    
    std::vector<std::complex<double>> orig(coefficients);
    std::vector<std::complex<double>> nule;
    int t(0);
    bool c(false);
    std::complex<double> x, w, v;
    
    for (int i=coefficients.size()-1; i>=1; i--) {
        t=1;
        c=false;
        while (!c && t<maxtrials) {
            x=RandomComplex(-10, 10, -10, 10);
            auto rez = Laguerre(coefficients, i, x, eps, maxiters);
            x=rez.first;
            c=rez.second;
            t++;
        }
        if (!c) throw std::logic_error("Convergence has not achieved");
        
        //Poliranje:
        auto rez = Laguerre(orig, orig.size()-1, x, eps, maxiters);
        if (rez.second) x=rez.first;
        
        if (std::fabs(x.imag())<=eps) x=x.real();
        nule.push_back(x);
        
        v=coefficients.at(i);
        for (int j=i-1; j>=0; j--) {
            w=coefficients.at(j);
            coefficients.at(j)=v;
            v=w+x*v;
        }
    }
    
    std::sort(nule.begin(), nule.end(), [](std::complex<double> c1, std::complex<double> c2){ return (c1.real()<c2.real() || (doubleequal(c1.real(), c2.real()) && c1.imag()<c2.imag())); });
    
    return nule;
}
std::vector<std::complex<double>>PolyRoots(std::vector<double> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10) {
    if (eps<=0 || maxiters<=0 || maxtrials<=0) throw std::domain_error("Invalid parameters");
    
    std::vector<double> orig(coefficients);
    std::vector<std::complex<double>> nule(coefficients.size());
    std::complex<double> x;
    int i(coefficients.size()-1), t(0);
    bool c(false);
    double v(0), a(0), b(0), u(0), w(0);
    
    while (i>=1) {
        t=1;
        c=false;
        while (!c && t<maxtrials) {
            x=RandomComplex(-10, 10, -10, 10);
            auto rez = Laguerre(coefficients, i, x, eps, maxiters);
            x=rez.first;
            c=rez.second;
            t++;
        }
        if (!c) throw std::logic_error("Convergence has not achieved");
        
        //Poliranje:
        auto rez = Laguerre(orig, orig.size()-1, x, eps, maxiters);
        if (rez.second) x=rez.first;
        
        if (std::fabs(x.imag())<=eps) {
            x=x.real();
            nule.at(i)=x;
            v=coefficients.at(i);
            for (int j=i-1; j>=0; j--) {
                w=coefficients.at(j);
                coefficients.at(j)=v;
                v=w+x.real()*v;
            }
            i--;
        }
        else {
            nule.at(i)=x;
            nule.at(i-1)=conj(x);
            a=2*x.real();
            b=std::abs(x)*std::abs(x);
            u=coefficients.at(i);
            v=coefficients.at(i-1)+a*u;
            for (int j=i-2; j>=0; j--) {
                w=coefficients.at(j);
                coefficients.at(j)=u;
                u=v;
                v=w+a*v-b*coefficients.at(j);
            }
            i=i-2;
        }
    }
    nule.erase(nule.begin());
    std::sort(nule.begin(), nule.end(), [](std::complex<double> c1, std::complex<double> c2){ return (c1.real()<c2.real() || (doubleequal(c1.real(), c2.real()) && c1.imag()<c2.imag())); });
    
    return nule;
}

template <typename FunTip>
void OgradiMinimum(FunTip f, double &x, double &h, double hmax, double lambda, double  &a, double &b) {
    
    while (std::fabs(h)<hmax) {
        if (f(x+h)<f(x)) {
            b=x+h;
            a=x-h;
        }
        else if (f(x-h)<f(x)) {
            b=x-h;
            a=b-h;
        }
        else {
            a=x-h;
            b=x+h;
            break;
        }
        x=b;
        h*=lambda;
    }
    if (std::fabs(h)>=hmax) throw std::logic_error("Minimum has not found");
}

template <typename FunTip>
double FindMinimum(FunTip f, double x0, double eps = 1e-8, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4) {
    if (eps<=0 || hinit<=0 || hmax<=0 || lambda<=0) throw std::domain_error("Invalid parameters");
    
    double h(hinit);
    double x(x0);
    double a(x-h), b(x+h);
    
    OgradiMinimum(f, x, h, hmax, lambda, a, b);
    
    double d, u, v, phi=(1+std::sqrt(5))/2, c(x);
    if (std::fabs(c-a)<std::fabs(b-c)) d=b-(b-c)/phi;
    else {
        d=c;
        c=a+(c-a)/phi;
    }
    u=f(c);
    v=f(d);
    
    while (std::fabs(b-a)>eps) {
        if (u<v) {
            b=d;
            d=c;
            c=a+(c-a)/phi;
            v=u;
            u=f(c);
        }
        else {
            a=c;
            c=d;
            d=b-(b-d)/phi;
            u=v;
            v=f(d);
        }
    }
    return (a+b)/2;
}

template <typename FunTip>
std::vector<std::pair<double, double>> RK4Integrator(FunTip f, double x0, double y0, double xmax, double h, double eps = 1e-8, bool adaptive = false) {
    std::vector<std::pair<double, double>> tacke;
    
    double x(x0), y(y0), step, u, v, w, delta;
    if (h<0) xmax=-xmax;
    if (!adaptive) tacke.push_back({x, y});
    
    while (std::fabs(x)<=xmax+eps) {
        
        if (!adaptive) {
            step=RK4Step(f, x, y, h);
            y=step;
            x=x+h;
            if (std::fabs(x)<=xmax+eps) tacke.push_back({x, y});
        }
        else {
            u=RK4Step(f, x, y, h/2);
            v=RK4Step(f, x+h/2, u, h/2);
            w=RK4Step(f, x, y, h);
            delta=std::fabs(w-v)/std::fabs(h);
            if (delta<=eps) {
                x=x+h;
                y=v;
                tacke.push_back({x, y});
            }
            h=h*std::min(5.0, 0.9*std::pow((eps/delta), 1./4));
        }
    }
    
    if(std::fabs(tacke.back().first)>xmax && adaptive) {
        if (h<0) xmax=-xmax;
        h=xmax-tacke[tacke.size()-1].first;
        u=RK4Step(f,x,y,h/2);
        v=RK4Step(f,x+h/2,u,h/2);
        tacke[tacke.size()-1]={xmax,v};
    }
    
    
    return tacke;
}

void FFT(std::vector<double> &x, std::vector<std::complex<double>> &xp, int N, int s=0, int d=0, int t=1) {
    if (N==1) xp[d]=x[s];
    else {
        FFT(x, xp, N/2, s, d, 2*t);
        FFT(x, xp, N/2, s+t, d+N/2, 2*t);
        std::complex<double> mi, w, u, v;
        mi={1, 0};
        w=std::exp(std::complex<double>(0, (-2*pi)/N));
        
        for(int k=d; k<=(d+N/2-1); k++) {
            u=xp[k];
            v=mi*xp[k+N/2];
            xp[k]=u+v;
            xp[k+N/2]=u-v;
            mi=mi*w;
        }
    }
}

void invFFT(std::vector<std::complex<double>> &xp, std::vector<std::complex<double>> &x, int N, int s=0, int d=0, int t=1) {
    if (N==1) x[d]=xp[s];
    else {
        invFFT(xp, x, N/2, s, d, 2*t);
        invFFT(xp, x, N/2, s+t, d+N/2, 2*t);
        std::complex<double> mi, w, u, v;
        mi={1, 0};
        w=std::exp(std::complex<double>(0, (2*pi)/N));
        
        for(int k=d; k<=(d+N/2-1); k++) {
            u=x[k];
            v=mi*x[k+N/2];
            x[k]=(u+v)/2.;
            x[k+N/2]=(u-v)/2.;
            
            mi=mi*w;
        }
    }
}

std::vector<double> LossyCompress(std::vector<double> data, int new_size) {
    if(new_size<=1 || new_size>data.size()) throw std::range_error("Bad new size");
    
    int test(data.size());
    while (test%2==0) test=test/2;
    if (test!=1)  throw std::range_error("Data size must be a power of two");
    
    std::vector<double> y(data.size());
    for (int i=0; i<data.size()/2; i++) y.at(i)=data.at(2*i);
    for (int i=data.size()/2; i<data.size(); i++) y.at(i)=data.at(2*(data.size()-i)-1);
    
    std::vector<std::complex<double>> ydft(y.size());
    FFT(y, ydft, y.size());
    
    
    std::vector<double> xdft(ydft.size());
    for(int i=0; i<ydft.size(); i++) {
        xdft[i]=( std::exp(std::complex<double>(0, (-i*pi)/(2*data.size()))) * ydft[i]).real();
    }
    
    xdft.resize(new_size);
    xdft[new_size-1]=data.size();
    return xdft;
}

std::vector<double> LossyDecompress(std::vector<double> compressed) {
    int N=compressed[compressed.size()-1];
    
    if (N<1 || N<compressed.size()) throw std::logic_error("Bad compressed sequence");
    
    int test(N);
    while (test%2==0) test=test/2;
    if (test!=1)  throw std::logic_error("Bad compressed sequence");
    
    
    compressed.at(compressed.size()-1)=0;
    for (int i=compressed.size(); i<N; i++) compressed.push_back(0);
    
    std::vector<std::complex<double>> y(N);
    std::vector<std::complex<double>> ytmp(N);
    ytmp[0]=compressed[0];
    for (int i=1; i<N; i++) {
        ytmp[i]=2.*std::exp(std::complex<double>(0, (pi*i)/(2*N)))*compressed[i];
    }
    
    invFFT(ytmp, y, N);
    
    for (int i=0; i<N; i++) {
        if (i%2==0) compressed[i]=y[i/2].real();
        else compressed[i]=y[N-(i+1)/2].real();
    }
    
    return compressed;
}


void testLossyCompressAndDecompress() {
    int brojac(0);
    
    std::vector<double> x1(32), y, x2;
    for (int i=0; i<32; i++) x1[i]=rand()%10;
    
    y=LossyCompress(x1, 32);
    x2=LossyDecompress(y);
    
    for (int i=0; i<32; i++) if (int(std::fabs(x1[i]-x2[i]))==0) brojac++;
    
    if (brojac==32) std::cout << "Test LossyCompressAndDecompress: OK" << std::endl;
    else std::cout << "Test LossyCompressAndDecompress: Greska" << std::endl;
}
void testRK4Integrator() {
    int brojac(0);
    
    auto dif=[](double x, double y){return 2*x+3*y-11;};
    //Pocetni uslov y(0)=1
    auto rj=[](double x){return 1./9*(-6*x-22*std::exp(3*x)+31);};
    
    //bez adaptacije:
    std::vector<std::pair<double, double>> v=RK4Integrator(dif, 0, 1, 2, 0.1);
    
    int n1(v.size());
    for (int i=0; i<n1; i++) {
        if (std::fabs(v[i].second-rj(v[i].first))<0.00001) brojac++;
    }
    
    
    //sa adaptacijom:
    v=RK4Integrator(dif, 0, 1, 2, 0.1, 1e-8, true);
    
    int n2(v.size());
    for (int i=0; i<n2; i++) {
        if (std::fabs(v[i].second-rj(v[i].first))<0.00001) brojac++;
    }
    
    if (brojac==n1+n2) std::cout << "Test RK4Integrator: OK" << std::endl;
    else std::cout << "Test RK4Integrator: Greska" << std::endl;
}
void testFindMinimum() {
    int brojac(0);
    auto f=[](double x){return x*x+1;};
    
    //Neispravni parametri:
    try {
        FindMinimum(f, 1, -1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    
    double min = FindMinimum(f, 1.2);
    if (std::fabs(min)<0.00001) brojac++;
    
    if (brojac==2) std::cout << "Test FindMinimum: OK" << std::endl;
    else std::cout << "Test FindMinimum: Greska" << std::endl;
}
void testPolyRoots() {
    int brojac(0);
    std::vector<std::complex<double>> v1({{169, 0}, {10, 0}, {1, 0}}); 
    std::vector<std::complex<double>> nule1({{-5, -12}, {-5, 12}});
    
    std::vector<double> v2({153, -145, -9, 1});
    std::vector<std::complex<double>> nule2({{-9, 0}, {1, 0}, {17, 0}});
    
    std::vector<std::complex<double>> v3({{-42, -9}, {29, -6}, {-8, 3}, 1});
    std::vector<std::complex<double>> nule3({{1, 2}, 3, {4, -5}});
    
    //Neispravni parametri:
    try {
        PolyRoots(v1, -1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    try {
        PolyRoots(v2, -1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    try {
        auto rez1 = PolyRoots(v1);
        auto rez2 = PolyRoots(v2);
        auto rez3 = PolyRoots(v3);
        
        if (isti(rez1, nule1)) brojac++; 
        if (isti(rez2, nule2)) brojac++;
        if (isti(rez3, nule3)) brojac++; 
    }
    catch(...) {}
    
    if (brojac==5) std::cout << "Test PolyRoots: OK" << std::endl;
    else std::cout << "Test PolyRoots: Greska" << std::endl;
}
void testNewtonRaphsonSolve() {
    int brojac(0);
    double rez(0);
    
    auto f1=[](double x){return log(x);};
    auto f1prim=[](double x){return 1./x;};
    
    auto f2=[](double x){return (x+2)*(x-7);};
    auto f2prim=[](double x){return x-5;};
    
    auto f3=[](double x){return std::sin(x);};
    auto f3prim=[](double x){return std::cos(x);};
    
    //Neispravni parametri:
    try {
        NewtonRaphsonSolve(f2, f2prim, 3, -10);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    
    try {
        rez=NewtonRaphsonSolve(f1, f1prim, 3, 1e-10, 0.5);
        if (std::fabs(rez-1)<0.00001) brojac++;
        
        rez=NewtonRaphsonSolve(f2, f2prim, 3, 1e-10, 0.5);
        if (std::fabs(rez+2)<0.00001 || std::fabs(rez-7)<0.00001) brojac++;
        
        rez=NewtonRaphsonSolve(f3, f3prim, 2.6);
        if (std::fabs(rez/pi)-1<0.00001) brojac++;
    }
    catch(...) {}
    
    
    if (brojac==4) std::cout << "Test NewtonRaphsonSolve: OK" << std::endl;
    else std::cout << "Test NewtonRaphsonSolve: Greska" << std::endl;
}
void testRiddersSolve() {
    int brojac(0);
    double a(0), b(1), rez(0);
    
    auto f1=[](double x){return log(x);};
    auto f2=[](double x){return (x+2)*(x-7);};
    auto f3=[](double x){return std::sin(x);};
    
    //Neispravni parametri:
    try {
        RiddersSolve(f2, a, b);
    }
    catch(std::range_error) { brojac++; }
    catch(...) {}
    
    try {
        a=-3; b=-1;
        RiddersSolve(f2, a, b, IllinoisSlavic, -1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {} 
    
    
    
    try {
        a=0.1; b=1.9;
        rez=RiddersSolve(f1, a, b);
        if (doubleequal(1, rez)) brojac++;
        
        a=-2.2; b=-1.8;
        rez=RiddersSolve(f2, a, b);
        if (std::fabs(rez+2)<0.00001) brojac++;
        
        
        a=1.3; b=4.99;
        rez=RiddersSolve(f3, a, b);
        if (std::fabs(rez-pi)<0.001) brojac++;
    }
    catch(...) {}
    
    
    if (brojac==5) std::cout << "Test RiddersSolve: OK" << std::endl;
    else std::cout << "Test RiddersSolve: Greska" << std::endl;
}
void testRegulaFalsiSolve() {
    int brojac(0);
    double a(0), b(1), rez(0);
    
    auto f1=[](double x){return log(x);};
    auto f2=[](double x){return (x+2)*(x-7);};
    auto f3=[](double x){return std::sin(x);};
    
    //Neispravni parametri:
    try {
        RegulaFalsiSolve(f2, a, b);
    }
    catch(std::range_error) { brojac++; }
    catch(...) {}
    
    try {
        a=-3; b=-1;
        RegulaFalsiSolve(f2, a, b, IllinoisSlavic, -1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {} 
    
    
    //Unmodified:
    try {
        a=-2.2; b=-1.8;
        rez=RegulaFalsiSolve(f2, a, b, Unmodified);
        if (std::fabs(rez+2)<0.0001) brojac++;
        
        a=0.001; b=18;
        rez=RegulaFalsiSolve(f1, a, b, Unmodified);
    }
    catch(std::logic_error) { brojac++; }
    catch(...) {}
    
    
    
    //Slavic:
    try {
        a=-2.2; b=-1.8;
        rez=RegulaFalsiSolve(f2, a, b, Slavic);
        if (doubleequal(-2, rez)) brojac++;
        
        
        a=0.2; b=5;
        rez=RegulaFalsiSolve(f1, a, b, Slavic);
        if (doubleequal(1, rez)) brojac++;
    }
    catch(...) {}
    
    
    //Illinois
    try {
        a=-2.2; b=-1.8;
        rez=RegulaFalsiSolve(f2, a, b, Illinois);
        if (std::fabs(rez+2)<0.0001) brojac++;
        
        
        a=2; b=4;
        rez=RegulaFalsiSolve(f3, a, b, Illinois);
        if (doubleequal(pi, rez)) brojac++;
        
        a=0.2; b=1.2;
        rez=RegulaFalsiSolve(f1, a, b, Illinois);
        if (doubleequal(1, rez)) brojac++;
    }
    catch(...) {}
    
    
    
    //IllinoisSlavic
    try {
        a=-2.2; b=-1.8;
        rez=RegulaFalsiSolve(f2, a, b, IllinoisSlavic);
        if (std::fabs(rez+2)<0.00001) brojac++;
        
        a=2; b=4;
        rez=RegulaFalsiSolve(f3, a, b, IllinoisSlavic);
        if (doubleequal(pi, rez)) brojac++;
        
        a=0.2; b=1.2;
        rez=RegulaFalsiSolve(f1, a, b, IllinoisSlavic);
        if (doubleequal(1, rez)) brojac++;
    }
    catch(...) {}
    
    
    if (brojac==12) std::cout << "Test RegulaFalsiSolve: OK" << std::endl;
    else std::cout << "Test RegulaFalsiSolve: Greska" << std::endl;
}
void testBracketRoot() {
    int brojac(0);
    bool rez;
    double a(0), b(0);
    
    auto f1=[](double x){return x*x+1;};
    auto f2=[](double x){return (x+2)*(x+2);};
    auto f3=[](double x){return x*x + 4*x - 45;};
    
    //Neispravni parametri:
    try {
        rez = BracketRoot(f1, 3, a, b, -1);
    }
    catch(std::domain_error) { brojac++; }
    catch(...) {}
    
    
    //Nema nulu:
    rez = BracketRoot(f1, 2, a, b);
    if (!rez) brojac++;
    
    //Visestruka nula:
    rez = BracketRoot(f2, -6, a, b);
    if (!rez) brojac++;
    
    //Dvije nule:
    rez = BracketRoot(f3, -1, a, b);
    if (rez) brojac++;
    if ((a<-9 && b>-9) || (a<5 && b>5)) brojac++;
    
    if (brojac==5) std::cout << "Test BracketRoot: OK" << std::endl;
    else std::cout << "Test BracketRoot: Greska" << std::endl;
}

int main ()
{
    testBracketRoot();
    testRegulaFalsiSolve();
    testRiddersSolve();
    testNewtonRaphsonSolve();
    testPolyRoots();
    testFindMinimum();
    testRK4Integrator();
    testLossyCompressAndDecompress();
	return 0;
}