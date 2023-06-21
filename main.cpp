#include <bits/stdc++.h>
#include "matrix.cpp"
#include <eigen3/Eigen/Dense>
using namespace std;

mt19937 rnd(chrono::steady_clock::now().time_since_epoch().count());
#define rd(n) uniform_int_distribution<int>(1,n)(rnd)

typedef std::complex<double> complx;
typedef std::vector<complx> complx_vector;
typedef matrix_lib::matrix m1;
typedef Eigen::MatrixXcd m2;
const int T=1000, mx=10;
const double eps=1e-4;

void gen_matrix(int n, int m, m1 &a, m2 &b){
    a.resize(n,m);
    complx qwq;
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++){
            qwq={2.0/mx*(rd(mx)-1),2.0/mx*(rd(mx)-1)};
            // qwq={1,0};
            a.a[i][j]=qwq;
            b(i,j)=qwq;
        }
}
bool isequal(complx x, complx y){
    // cerr<<x<<" "<<y<<" isequal\n";
    return fabs(x.real()-y.real())<eps && fabs(x.imag()-y.imag())<eps;
}
bool check(const m1 &a, const m2 &b){
    for(int i=0; i<a.n; i++)
        for(int j=0; j<a.m; j++)
            if(!isequal(a.a[i][j], b(i,j)))return 0;
    return 1;
}

int n,m,k;
m1 a,b,c,d;

void add_sub_test(){
    puts("begin add_sub_test...");
    for(int TT=1; TT<=T; TT++){
        n=rd(10); m=rd(10);
        m2 x(n,m),y(n,m),z(n,m);
        gen_matrix(n,m,a,x);
        gen_matrix(n,m,b,y);
        c=a+b; z=x+y;
        // c.show(); cout<<z<<"\n"; cout<<isequal(c.a[0][0],z(0,0))<<" "<<check(c,z)<<"\n";
        assert(check(c,z));
        c=a-b; z=x-y;
        assert(check(c,z));
    }
    puts("finish add_sub_test...\n");
}

void mul_test(){
    puts("begin mul_test...");
    for(int TT=1; TT<=T; TT++){
        n=rd(10);m=rd(10);k=rd(10);
        m2 x(n,m), y(m,k), z(n,k);
        gen_matrix(n,m,a,x);
        gen_matrix(m,k,b,y);
        c=a*b; z=x*y;
        assert(check(c,z));
    }
    puts("finish mul_test...\n");
}

void del_test(){
    puts("begin det_test...");
    for(int TT=1; TT<=T; TT++){
        n=rd(10);
        m2 x(n,n);
        gen_matrix(n,n,a,x);
        // cerr<<a.det()<<" "<<x.determinant()<<" del\n";
        assert(isequal(a.det(),x.determinant()));
    }
    puts("finish det_test...\n");
}

bool comp_complx(const complx &x, const complx &y){
    return matrix_lib::sgn(x.real()-y.real())?x.real()<y.real():x.imag()<y.imag();
}

void eigenvalue_test(){
    puts("begin eigenvalue_test...");
    complx_vector qwq1;
    Eigen::VectorXcd qwq2;
    for(int TT=1; TT<=T; TT++){
        if((TT-1)%10==0) cerr<<"Running on Eigenvalue testcase "<<TT<<" ...\n";
        n=rd(10);
        // n=4;
        m2 x(n,n);
        gen_matrix(n,n,a,x);
        qwq1=a.eigenvalue();
        // Eigen::EigenSolver<Eigen::MatrixXcd> es(x);
        qwq2=x.eigenvalues();
        sort(qwq1.begin(),qwq1.end(),comp_complx);
        sort(qwq2.data(), qwq2.data()+qwq2.size(), comp_complx);
        for(int i=0; i<n; i++){
            if(!isequal(qwq1[i],qwq2[i])){
                for(auto v:qwq1) cerr<<v<<" ";
                cerr<<" mine\n";
                for(auto v:qwq2) cerr<<v<<" ";
                cerr<<" eigen\n";
            }
            assert(isequal(qwq1[i], qwq2[i]));
        }
    }
    puts("finish eigenvalue_test...\n");
}

void SVD_test(){
    puts("begin SVD_test...");
    for(int TT=1; TT<=T; TT++){
        if(TT%10==0) cerr<<"Running on SVD testcase "<<TT<<" ...\n";
        // n=2; m=n;
        n=rd(10); m=rd(10);
        m2 x(n,m);
        gen_matrix(n,m,a,x);

        // Eigen::JacobiSVD<m2> svd(x, Eigen::ComputeThinU | Eigen::ComputeThinV);
        // cout<<svd.matrixU()<<"\n";
        // cout<<svd.matrixV()<<"\n";
        // cout<<svd.singularValues()<<"\n\n";

        a.SVD_decomp(b,c,d);
        b=b.transpose();
        for(int i=0; i<n; i++)
            for(int j=i+1; j<n; j++){
                if(!isequal(matrix_lib::vec_dot(b.a[i],b.a[j]), complx(0,0))){
                    cerr<<matrix_lib::vec_dot(b.a[i],b.a[j])<<"\n";
                }
                assert(isequal(matrix_lib::vec_dot(b.a[i],b.a[j]), complx(0,0)));
            }
        for(int i=0; i<m; i++)
            for(int j=i+1; j<m; j++)
                assert(isequal(matrix_lib::vec_dot(d.a[i],d.a[j]), complx(0,0)));
    }
    puts("finish SVD_test...");
}

int main(){
    // add_sub_test();
    // mul_test();
    // del_test();
    eigenvalue_test();
    // SVD_test();
    return (0-0);
}