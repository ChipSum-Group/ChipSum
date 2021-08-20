


##快速开始

### 环境准备
- 使用Linux或macOS操作系统.

- `ChipSum`依托`kokkos`和`kokkos-kernels`编译, 编译需要`Cmake>3.16`，推荐使用`3.20.0`。
```
# get ChipSum
git clone https://github.com/CDCS-HPC/ChipSum.git

# get kokkos and kokkos-kernels
cd ChipSum/tpl
git clone https://github.com/kokkos/kokkos.git
git clone https://github.com/kokkos/kokkos-kernels.git
```

### 编译
- `ChipSum`使用了一个python脚本帮助完成编译过程。笔者分别在`AMD VEGA906`和`NVIDIA 2080ti` 架构上完成编译。
```
# AMD Vega906/900  
python3 setup.py arch=VEGA906 compiler=/path/to/your/rocm(i.e. 4.0.1)/bin/hipcc hip=/path/to/your/rocm(i.e. 4.0.1)

# NVIDIA 2080ti
python3 setup.py arch=Turing75 cuda=/Path/To/Your/Cuda

# for Volta72
python3 setup.py cuda=/Path/To/Your/Cuda arch=Volta72

# for MX350（Pascal61）
python3 setup.py cuda=/Path/To/Your/Cuda arch=Pascal61
```
- 若想指定安装目录，可以使用`prefix`参数
```
mkdir /anywhere/you/like
python3 setup.py prefix=/anywhere/you/like
```
- 若想指定核编译，可以使用`j`参数
```
python3 setup.py j=32
```
- 第一次编译时需编译kokkos和kokkos-kernels，耗时较久。后续使用时仅编译ChipSum内容，耗时很快。

### 验证安装

 `ChipSum`在`test.cpp`中提供了一个简单的用例示范，默认路径编译完成后，可以在`./build`中查看编译结果。
```
# default path
cd ./build
./ChipSum
```
预期输出：
```
Kokkos::OpenMP::initialize WARNING: OMP_PROC_BIND environment variable not set
  In general, for best performance with OpenMP 4.0 or better set OMP_PROC_BIND=spread and OMP_PLACES=threads
  For best performance with OpenMP 3.1 set OMP_PROC_BIND=true
  For unit testing set OMP_PROC_BIND=false
vector_0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
3.3
scalar_1: [0, 3.3, 6.6, 9.9, 13.2, 16.5, 19.8, 23.1, 26.4, 29.7]
```