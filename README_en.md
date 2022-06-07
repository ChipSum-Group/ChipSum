<img width="578" alt="chipsum-logo" src="https://user-images.githubusercontent.com/3610126/123236443-61fb0a80-d50f-11eb-8473-0c158979b9f7.png">

# Introduction to ChipSum
- ChipSum is a framework intigrating HPC libraries using C++ template programming. 
- Chipsum is oriented to develop a performance potability platform for diverse heterogeneous computing devices, including NVIDIA/AMD's manycore GPU and Intel's multicore CPU.
- ChipSum's core code is a wrap of vector/matrix operations based on 3rd party libraries.
- ChipSum's core code can derive many child codes for special issues, currently including foundation codes ChipSum.Solver, ChipSum.AI, and application code ChipSum.MD ...

# Architecture of ChipSum
#### ChipSum's ecosystem
|Level Defination|Description|Contents|
|:--:|:--:|:--:|
|Level 3|ChipSum's application code|[ChipSum.MD](https://github.com/CDCS-HPC/ChipSum.MD)|
|Level 2|ChipSum's foundation code|[ChipSum.Solver](https://github.com/CDCS-HPC/ChipSum.Solver), [ChipSum.AI](https://github.com/CDCS-HPC/ChipSum.AI)|
|Level 1|ChipSum's core code|[ChipSum](https://github.com/CDCS-HPC/ChipSum)|

ChipSum is a main framework intigrating HPC libraries using C++ template programming.


# Quick Start

[Chinese](./README_CN.md)

### Environment

- Available on `Linux` or `macOS`. 

`ChipSum` compile with  `kokkos` and `kokkos-kernels`. `Cmake` > `3.16`, `3.20.0` is recommended.

```
# get ChipSum
git clone https://github.com/CDCS-HPC/ChipSum.git

# get kokkos and kokkos-kernels
cd ChipSum/tpl
git clone https://github.com/kokkos/kokkos.git
git clone https://github.com/kokkos/kokkos-kernels.git
```


### Compile

- `ChipSum` can be compiled automatically by a python3 script `setup.py`. We finished compiling on both `AMD Vega 906` and `NVIDIA 2080TI`

```
# for AMD Vega906/900  
python3 setup.py arch=VEGA906 compiler=/path/to/your/rocm(i.e. 4.0.1)/bin/hipcc hip=/path/to/your/rocm(i.e. 4.0.1)

# for NVIDIA 2080ti
python3 setup.py arch=Turing75 cuda=/Path/To/Your/Cuda

# for Volta72
python3 setup.py cuda=/Path/To/Your/Cuda arch=Volta72

# for MX350(Pascal61)
python3 setup.py cuda=/Path/To/Your/Cuda arch=Pascal61
```

- If you want to specify an install path, use `prefix`
```
mkdir /any/path/you/like
python3 setup.py prefix=/any/path/you/like 
```

- If you want to specify compile cores, use `j`

```
python3 setup.py j=32 
```

PS: it takes a long time when you compile with kokkos and kokkos-kernels for the first time and it will be fast when you compile ChipSum application after that.


### Verify Installation
`ChipSum` provides a toy example in `test.cpp`. You can find the application in `./build` if compiled by default.

```
cd ./build
./ChipSum
```
Expecting output:

```
Kokkos::OpenMP::initialize WARNING: OMP_PROC_BIND environment variable not set
  In general, for best performance with OpenMP 4.0 or better set OMP_PROC_BIND=spread and OMP_PLACES=threads
  For best performance with OpenMP 3.1 set OMP_PROC_BIND=true
  For unit testing set OMP_PROC_BIND=false

(test.cpp output)

```

