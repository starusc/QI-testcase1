#include <bits/stdc++.h>

namespace matrix_lib{

#define O(x) std::cerr<<(#x)<<" : "<<x<<"\n";
typedef std::complex<double> complx;
typedef std::vector<complx> complx_vector;
typedef std::vector<complx_vector> complx_matrix;
const double eps=1e-8;
const int max_iter=1000;

class matrix{
    public:
        int n, m;
        complx_matrix a;
    public:
        matrix(){}
        matrix(int nn, int mm, const complx_vector &b);
        void resize(int nn, int mm, int fl=0);
        void init(int nn, int mm, const complx_vector &b);
        void show();
        bool operator ==(const matrix &b);
        matrix operator +(const matrix &b);
        matrix operator -(const matrix &b);
        matrix operator *(const matrix &b);
        complx_vector operator *(const complx_vector &b);
        matrix transpose();
        matrix conjugate_transpose();
        void QR_decomp(matrix &q, matrix &r);
        complx det();
        complx_vector eigenvalue();
        void SVD_decomp(matrix &u, matrix &d, matrix &v);
};

matrix::matrix(int nn, int mm, const complx_vector &b){
    matrix::init(nn, mm, b);
}

void matrix::resize(int nn, int mm, int fl){
    n=nn; m=mm;
    a.resize(n);
    for(int i=0; i<n; i++) a[i].resize(m);
    if(fl){
        for(int i=0; i<n; i++)
            for(int j=0; j<m; j++) a[i][j]=complx(0,0);
    }
}

void matrix::init(int nn, int mm, const complx_vector &b){
    n=nn; m=mm;
    resize(n,m);
    for(int i=0, x=0; i<n; i++)
        for(int j=0; j<m; j++) a[i][j]=b[x++];
}

void matrix::show(){
    printf("%d row %d columns:\n", n, m);
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            printf("%.6lf+%.6lfi ", a[i][j].real(), a[i][j].imag());
        }
        printf("\n");
    }
    printf("showing matrix finished.\n");
}

int sgn(double x){return fabs(x)<eps?0:(x>0?1:-1);}
int sgn(complx x){return sgn(fabs(x.real())+fabs(x.imag()));}

bool matrix::operator ==(const matrix &b){
    if(n!=b.n || m!=b.m) return 0;
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++)
            if(sgn(a[i][j]-b.a[i][j])) return 0;
    return 1;
}

matrix matrix::operator +(const matrix &b){
    if(n!=b.n || m!=b.m){
        std::cerr<<"Error! Add operation: dimension do not match.\n";
        exit(0);
    }
    matrix c; c.resize(n,m);
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++)
            c.a[i][j]=a[i][j]+b.a[i][j];
    return c;
}

matrix matrix::operator -(const matrix &b){
    if(n!=b.n || m!=b.m){
        std::cerr<<"Error! Sub operation: dimension do not match.\n";
        exit(0);
    }
    matrix c; c.resize(n,m);
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++)
            c.a[i][j]=a[i][j]-b.a[i][j];
    return c;
}

matrix matrix::operator *(const matrix &b){
    if(m!=b.n){
        std::cerr<<"Error! Mul operation: dimension do not match.\n";
        exit(0);
    }
    matrix c; c.resize(n,b.m);
    for(int i=0; i<n; i++)
        for(int j=0; j<b.m; j++){
            c.a[i][j]=0;
            for(int k=0; k<m; k++)
                c.a[i][j]+=a[i][k]*b.a[k][j];
        }
    return c;
}

complx_vector matrix::operator *(const complx_vector &b){
    if(m!=(signed)b.size()){
        std::cerr<<"Error! Vector mul operation: dimension do not match.\n";
        exit(0);
    }
    complx_vector c; complx x;
    for(int i=0; i<n; i++){
        x=complx(0,0);
        for(int j=0; j<m; j++) x+=a[i][j]*b[j];
        c.push_back(x);
    }
    return c;
}

matrix matrix::transpose(){
    matrix c; c.resize(m,n);
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++)
            c.a[j][i]=a[i][j];
    return c;
}

matrix matrix::conjugate_transpose(){
    matrix c; c.resize(m,n);
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++)
            c.a[j][i]=conj(a[i][j]);
    return c;
}

complx matrix::det(){
    if(n!=m){
        std::cerr<<"Error! Det operation: matrix is not square.\n";
        exit(0);
    }
    complx ans(1,0), tmp;
    complx_matrix d(a);
    for(int i=0; i<n; i++){
        if(!sgn(d[i][i])){
            for(int j=i+1; j<n; j++) if(sgn(d[j][i])){
                swap(d[i], d[j]);
                ans=-ans;
                break;
            }
        }
        if(!sgn(d[i][i])){
            ans={0,0};
            break;
        }
        ans*=d[i][i];
        for(int j=i+1; j<n; j++){
            tmp=d[j][i]/d[i][i];
            for(int k=i; k<n; k++){
                d[j][k]-=tmp*d[i][k];
            }
        }
    }
    return ans;
}

complx vec_dot(const complx_vector &x, const complx_vector &y){
    complx ret(0,0);
    for(int len=x.size(), i=0; i<len; i++) ret+=x[i]*conj(y[i]);
    return ret;
}
double l2norm(const complx_vector &x){return sqrt(vec_dot(x,x).real());}
void vec_add(complx_vector &x, const complx_vector &v){
    for(int len=x.size(), i=0; i<len; i++) x[i]+=v[i];
}
void vec_mul(complx_vector &x, complx v){
    for(int len=x.size(), i=0; i<len; i++) x[i]*=v;
}
void vec_out(const complx_vector &x){
    for(auto v:x) std::cerr<<v<<" ";
    std::cerr<<"vecout\n";
}

void matrix::QR_decomp(matrix &q, matrix &r){
    complx_vector tmp;
    q=(*this).transpose(); r.resize(n,n);
    // q.show();
    for(int i=0; i<m; i++){
        r.a[i][i]=l2norm(q.a[i]);
        vec_mul(q.a[i], complx(1,0)/r.a[i][i]);
        for(int j=i+1; j<m; j++){
            r.a[i][j]=vec_dot(q.a[i], q.a[j]);
            tmp=q.a[i]; vec_mul(tmp, -r.a[i][j]);
            vec_add(q.a[j], tmp);
        }
        for(int j=0; j<i; j++) r.a[i][j]=0;
    }
    q=q.transpose();
    // q.show(); r.show(); (q*r).show();
    //check
    assert(q*r==(*this));
}

complx_vector matrix::eigenvalue(){
    if(n!=m){
        std::cerr<<"Error! Eigenvalue operation: matrix is not square.\n";
        exit(0);
    }
    complx tmp(0,0);
    matrix ak(*this), q, r;
    for(int T=1; T<=max_iter; T++){
        // tmp=ak.a[n-1][n-1];
        // for(int i=0; i<n; i++) ak.a[i][i]-=tmp;
        ak.QR_decomp(q, r);
        ak=r*q;
        // for(int i=0; i<n; i++) ak.a[i][i]+=tmp;
    }
    complx_vector ans;
    for(int i=0; i<n; i++) ans.push_back(ak.a[i][i]);
    return ans;
}

matrix rotate(int n, int p, int q, double theta){
    matrix c; c.resize(n,n,1);
    for(int i=0; i<n; i++) c.a[i][i]=complx(1,0);
    c.a[p][p]=c.a[q][p]=complx(cos(theta),-sin(theta))/sqrt(2.0);
    c.a[p][q]=complx(-cos(theta),-sin(theta))/sqrt(2.0);
    c.a[q][q]=complx(cos(theta),sin(theta))/sqrt(2.0);
    return c;
}

void matrix::SVD_decomp(matrix &u, matrix &d, matrix &v){
    matrix b(*this), R, J;
    double phi1, phi2;
    // todo: check n m
    b=b.conjugate_transpose()*b; // n>=m b:m*m
    u.resize(0,n); d.resize(n,m,1); v.resize(m,m); J.resize(m,m,1);
    for(int i=0; i<m; i++) J.a[i][i]=complx(1,0);

    for(int T=1, cnt; T<=max_iter; T++){
        cnt=0;
        for(int i=0; i<m; i++)
            for(int j=i+1; j<m; j++)if(sgn(b.a[i][j])){
                phi1=atan2(b.a[i][j].imag(), b.a[i][j].real());
                phi2=atan2(norm(b.a[i][j])*2, (b.a[i][i]-b.a[j][j]).real());
                R=rotate(m,i,j,phi2/2)*rotate(m,i,j,phi1/2-M_PI/4);
                b.show();
                b=R*b*R.conjugate_transpose();
                b.show(); R.show(); std::cerr<<"\n";
                //check
                J=J*R;
                ++cnt;
            }
        if(!cnt)break;
    }
    b.show(); J.show();
    
    std::vector<int> c; c.resize(m);
    for(int i=0; i<m; i++) c[i]=i;
    std::sort(c.begin(), c.end(), [&b](int x, int y){
        return b.a[x][x].real()>b.a[y][y].real();
    });
    for(int i=0 ,x; i<m; i++){
        x=c[i];
        d.a[i][i]=sqrt(b.a[x][x].real());
        v.a[i]=J.a[x];
        if(sgn(d.a[i][i])){
            vec_out(v.a[i]);vec_out((*this)*v.a[i]);
            u.a.push_back((*this)*v.a[i]); ++u.n;
            // for(int j=0; j<n; j++){
            //     std::cerr<<((*this)*v.a[i])[j]<<" "<<u.a[i][j]<<"\n";
            // }
            std::cerr<<l2norm(u.a[i])<<" "<<d.a[i][i]<<" ";
            vec_mul(u.a[i], complx(1,0)/d.a[i][i]);
            std::cerr<<l2norm(u.a[i])<<"\n";
        }
    }
    //expand U
    //check
    for(int i=0; i<u.n; i++){
        assert(sgn(l2norm(u.a[i])-1)==0);
        for(int j=i+1; j<u.n; j++)assert(sgn(vec_dot(u.a[i],u.a[j]))==0);
    }
    complx_vector tmp,qwq; tmp.resize(n);
    for(int i=0; i<n; i++){
        if(u.n==n)break;
        for(int j=0; j<n; j++) tmp[j]=complx(0,0);
        tmp[i]=complx(1,0);
        for(int j=0; j<u.n; j++){
            qwq=u.a[j];
            vec_mul(qwq, -vec_dot(tmp,u.a[j]));
            vec_add(tmp, qwq);
            // for(auto vv:tmp) std::cerr<<vv<<" ";
            // std::cerr<<"\n";
        }
        if(sgn(l2norm(tmp))){
            vec_mul(tmp,1/l2norm(tmp));
            u.a.push_back(tmp); ++u.n;
        }
    }
    if(u.n!=n) {u.show(); O(u.n);}
    assert(u.n==n);
    u=u.transpose(); // transpose忘记传给u了
    v=v.transpose().conjugate_transpose();

    //check
    // (*this).show(); u.show(); d.show(); v.show(); (u*d*v).show();
    assert(u*d*v==(*this));
}

}