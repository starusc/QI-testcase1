#include <bits/stdc++.h>

namespace matrix_lib{

#define O(x) std::cerr<<(#x)<<" : "<<x<<"\n";
typedef std::complex<double> complx;
typedef std::vector<complx> complx_vector;
typedef std::vector<complx_vector> complx_matrix;
const double eps=1e-8;
const int max_iter=10000;

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
        matrix operator *(const complx &x);
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
    if(fl&1){
        for(int i=0; i<n; i++)
            for(int j=0; j<m; j++) a[i][j]=complx(0,0);
    }
    if(fl&2){
        for(int i=0; i<n; i++) a[i][i]=complx(1,0);
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

matrix matrix::operator *(const complx &x){
    matrix c(*this);
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++) c.a[i][j]*=x;
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
void vec_clear(complx_vector &x, int n){
    x.resize(n);
    for(int i=0; i<n; i++) x[i]=complx(0,0);
}
bool vec_normalize(complx_vector &x){
    double qwq=l2norm(x);
    if(sgn(qwq)) {vec_mul(x,1/qwq); return true;}
    else {vec_clear(x,x.size()); return false;}
}

bool householder(const matrix &x, matrix &h, int i, int k){
    complx_vector w; vec_clear(w,x.n);
    h.resize(x.n,x.n,3);
    for(int j=k; j<x.n; j++) w[j]=x.a[j][i];
    if(sgn(norm(w[k]))) w[k]-=l2norm(w)/sqrt(norm(w[k]))*w[k]; // 复数与实数不一样
    else w[k]-=l2norm(w);
    if(!vec_normalize(w))return false;
    for(int u=i; u<x.n; u++)
        for(int v=i; v<x.n; v++) h.a[u][v]-=complx(2,0)*w[u]*conj(w[v]);
    return true;
}

void hessenberg(matrix &x){
    matrix h;
    for(int i=0; i<x.n-2; i++){
        if(!householder(x,h,i,i+1)) continue;
        x=h*x*h.conjugate_transpose();
    }
    //check
    // for(int i=0; i<x.n; i++)
    //     for(int j=0; j<i-1; j++)assert(sgn(x.a[i][j])==0);
    // x.show();
}

void matrix::QR_decomp(matrix &q, matrix &r){ // householder
    matrix h;
    r=*this; q.resize(n,n,3);
    for(int i=0; i<n-1; i++){
        if(!householder(r,h,i,i)) continue;
        r=h*r;
        q=q*h.conjugate_transpose();
    }

    //check
    // for(int i=0; i<n; i++){
    //     assert(sgn(l2norm(q.a[i])-1)==0);
    //     for(int j=i+1; j<i; j++)
    //         assert(sgn(vec_dot(q.a[i],q.a[j]))==0);
    // }
    // if(!(q*r==(*this))) {std::cerr<<"\n"; q.show(); r.show(); (q*r).show(); (*this).show();}
    assert(q*r==(*this));
}

complx_vector matrix::eigenvalue(){
    if(n!=m){
        std::cerr<<"Error! Eigenvalue operation: matrix is not square.\n";
        exit(0);
    }
    matrix ak(*this), M, q, r;
    complx_vector ans;
    complx tr, dt;
    int nn=n;
    hessenberg(ak);
    // ak.show();
    for(int T=1; T<=max_iter; T++){
        if(n==0)break;
        if(n==1){ans.push_back(ak.a[0][0]);n=0;break;}
        // M=H^2-(mu1+mu2)H+mu1mu2I mu1+mu2=-b/a mu1mu2=c/a
        /*tr=ak.a[n-2][n-2]+ak.a[n-1][n-1];
        dt=ak.a[n-2][n-2]*ak.a[n-1][n-1]-ak.a[n-2][n-1]*ak.a[n-1][n-2];
        M=ak*ak-ak*tr;
        for(int i=0; i<n; i++)M.a[i][i]+=dt;
        M.QR_decomp(q, r);
        ak=q.conjugate_transpose()*ak*q;
        if(sgn(q.a[n-1][n-1]-complx(1,0))==0&&sgn(q.a[n-2][n-2]-complx(1,0))==0){
            ans.push_back((tr+sqrt(tr*tr-complx(4,0)*dt))/complx(2,0));
            ans.push_back((tr-sqrt(tr*tr-complx(4,0)*dt))/complx(2,0));
            n-=2;
            ak.resize(n,n);
        }
        while(n>1&&sgn(ak.a[n-1][n-2])==0){
            ans.push_back(ak.a[n-1][n-1]);
            --n;
            ak.resize(n,n);
        }*/
        // q.show(); ak.show(); std::cerr<<"\n";
        ak.QR_decomp(q,r); ak=r*q;
    }
    // O(n); ak.show();
    for(int i=0; i<n; i++) ans.push_back(ak.a[i][i]);
    n=nn;
    return ans;
}

matrix rotate(int n, int p, int q, double y, double x){ // tan(2theta)=y/x
    x/=sqrt(y*y+x*x); // cos(2theta)
    double s=sqrt((1-x)/2), c=sqrt((1+x)/2);
    if(sgn(y)<0) s=-s;
    
    matrix ret; ret.resize(n,n,3);
    ret.a[p][p]=ret.a[q][p]=complx(c,-s)/sqrt(2.0);
    ret.a[p][q]=complx(-c,-s)/sqrt(2.0);
    ret.a[q][q]=complx(c,s)/sqrt(2.0);
    return ret;
}

void matrix::SVD_decomp(matrix &u, matrix &d, matrix &v){
    matrix A(*this), b, R, J; int fl=0;
    if(n<m){
        fl=1;
        std::swap(n,m);
        A=A.conjugate_transpose();
    }

    b=A.conjugate_transpose()*A; // n>=m b:m*m
    // b.show();
    u.resize(0,n); d.resize(n,m,1); v.resize(m,m); J.resize(m,m,3);

    for(int T=1, cnt; T<=max_iter; T++){
        cnt=0;
        for(int i=0; i<m; i++)
            for(int j=i+1; j<m; j++)if(sgn(b.a[i][j])){
                R=rotate(m,i,j,sqrt(norm(b.a[i][j]))*2, (b.a[i][i]-b.a[j][j]).real())
                    *rotate(m,i,j,-b.a[i][j].real(), b.a[i][j].imag()); // norm 没开根
                b=R*b*R.conjugate_transpose();
                // b.show(); R.show(); std::cerr<<"\n";
                J=J*R.conjugate_transpose();
                ++cnt;
            }
        if(!cnt)break;
    }
    // b.show(); J.show(); (J*b*J.conjugate_transpose()).show();
    J=J.transpose();
    
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
            u.a.push_back(A*v.a[i]); ++u.n;
            // std::cerr<<l2norm(u.a[i])<<" "<<d.a[i][i]<<" ";
            vec_mul(u.a[i], complx(1,0)/d.a[i][i]);
        }
    }

    //expand U
    complx_vector tmp,qwq; tmp.resize(n);
    for(int i=0; i<n; i++){
        if(u.n==n)break;
        for(int j=0; j<n; j++) tmp[j]=complx(0,0);
        tmp[i]=complx(1,0);
        for(int j=0; j<u.n; j++){
            qwq=u.a[j];
            vec_mul(qwq, -vec_dot(tmp,u.a[j]));
            vec_add(tmp, qwq);
        }
        if(vec_normalize(tmp)){
            u.a.push_back(tmp); ++u.n;
        }
    }

    assert(u.n==n);
    u=u.transpose(); // transpose忘记传给u了
    v=v.transpose().conjugate_transpose();
    if(fl){
        std::swap(n,m);
        J=v.conjugate_transpose();
        v=u.conjugate_transpose();
        u=J;
        d=d.transpose();
    }
    // (*this).show(); u.show(); d.show(); v.show(); (u*d*v).show();
    assert(u*d*v==(*this));
}

}