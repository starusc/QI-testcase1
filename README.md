## 运行说明文档

这是一个复数矩阵基本运算的函数集，用 C++ 语言编写，支持矩阵加法，矩阵乘法（乘矩阵、向量、常数），求行列式，求特征值，求SVD分解。

项目中 `matrix.cpp` 是函数集，`main.cpp` 是用于测试的程序（采用随机造数据，然后与 `Eigen` 库进行比对来测试）

### 初始化一个矩阵

`matrix qwq` 可以定义一个矩阵。

初始化传入三个参数 `int n, int m, vector<complex<double>> b` 分别表示矩阵的行数，列数，和矩阵的值，其中矩阵第 $i$ 行 $j$ 列（从 $0$ 开始标号）的值存储在 $b[i*m+j]$ 中。可以使用析构函数或函数 `void init(int nn, int mm, const complx_vector &b)` 初始化矩阵。

`void resize(int nn, int mm)` 可以重定义矩阵的大小。 

### 矩阵运算

判断矩阵相等，计算矩阵加法，计算矩阵减法，计算矩阵乘法（可以右乘矩阵、向量、常数）功能分别使用运算符 `==, +, -, *`。

`matrix transpose()` 返回当前矩阵的转置。`matrix conjugate_transpose()` 返回当前矩阵的共轭转置。`complx det()` 返回当前矩阵的行列式。

`void SVD_decomp(matrix &u, matrix &d, matrix &v)` 向此函数传入参数 $u,d,v$，函数返回时 $u,d,v$ 将被修改为当前矩阵的 SVD 分解，即当前矩阵 $=udv$。该函数采用双边 `Jacobi` 旋转迭代算法，参考 [wikipedia1](https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm), [wikipedia2](https://en.wikipedia.org/wiki/Jacobi_method_for_complex_Hermitian_matrices)。

`void QR_decomp(matrix &q, matrix &r)` 向该函数传入参数 $q,r$，函数返回时 $q,r$ 将被修改为当前矩阵的 QR 分解，即当前矩阵 $=qr$。该函数采用 `householder` 算法。

`complx_vector eigenvalue()` 返回当前矩阵的所有特征值。总体采用 [QR 算法](https://en.wikipedia.org/wiki/QR_algorithm)，添加两个优化。第一个一开始将初始矩阵转化成上海森堡矩阵。第二个为 [`double shift`](https://web.stanford.edu/class/cme335/lecture5) 优化，但此优化并未产生理想中加快收敛的效果，没有放到最终版本。

### 总结

这次实现的都是比较 general 的算法，加速空间还很大。第一可以优化代码，比如稀疏矩阵的乘法单独写可以快很多；第二可以针对不同的情况采用不同的算法，比如矩阵的大小，疏密，是否为对称矩阵等。这次实现的函数除了最后一个 `eigenvalue()` 外别的都能顺利通过小数据的检验；`eigenvalue()` 会在极少数情况上无法在设定的迭代次数内收敛，导致误差比较大。

想说的话：我第一眼看到这个项目的时候觉得很简单，心想这么常用的东西网上应该遍地都是。但是到了真正开始做的时候才发现资源比我想象中少很多，特别是不像学习的时候一样，不会写可以学习别人的代码然后掌握。刚开始的时候看到网上有很多很多算法，不知道该写啥。后来静下心来评估了各个算法的思路之后，选择了我认为性价比最高的。没找到啥复数的资料，所以写的过程中被实数变复数要做出的改变坑了好多次，一些看似正确的东西在出错了之后，仔细一推才发现转到复数域要做出改变。写的过程还挺有趣的，虽然自己菜菜的，但还是会成长起来的。