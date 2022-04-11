# 快速开始
## 一、编译环境准备
1. Linux系统

目前`ChipSum`仅支持Linux系统

2. CMake 

 `ChipSum`以及第三方库`kokkos`和`kokkos-kernels`的编译, 都需要`Cmake>3.16`，推荐使用`3.20.0`。

3. 获取`ChipSum`代码

```
# get ChipSum
git clone https://gitee.com/chip-sum/ChipSum.git

# get kokkos and kokkos-kernels
cd ChipSum/

git submodule init
git submodule update
```
## 二、编译
为了方便用户编译，`ChipSum`使用了一个python脚本帮助完成编译过程。已分别在`AMD VEGA906`和`NVIDIA 2080ti`等架构上完成编译并运行。

1. `AMD VEGA906`

```
# AMD Vega906/900  
python3 setup.py arch=VEGA906 compiler=/path/to/your/rocm(i.e. 4.0.1)/bin/hipcc hip=/path/to/your/rocm(i.e. 4.0.1)
```
2. `NVIDIA 2080ti`

```
# NVIDIA 2080ti
python3 setup.py arch=Turing75 cuda=/Path/To/Your/Cuda
```
3. `Volta72`
```
# for Volta72
python3 setup.py cuda=/Path/To/Your/Cuda arch=Volta72
```
4. `MX350(Pascal61)`
```
# for MX350（Pascal61）
python3 setup.py cuda=/Path/To/Your/Cuda arch=Pascal61
```
5. 其他

若想指定安装目录，可以使用`prefix`参数
```
mkdir /anywhere/you/like
python3 setup.py prefix=/anywhere/you/like
```
 若想指定核编译，可以使用`j`参数
```
python3 setup.py j=32
```
注意：第一次编译时需编译kokkos和kokkos-kernels，耗时较久。后续使用时仅编译ChipSum内容，耗时很快。

## 三、算例

 `ChipSum`在`test.cpp`中提供了一个简单的用例示范，默认路径编译完成后，可以在`./build`中查看编译结果。
```
# default path
cd ./build
./ChipSum
```
预期输出：
```
vector_0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
3.3
scalar_1: [0, 3.3, 6.6, 9.9, 13.2, 16.5, 19.8, 23.1, 26.4, 29.7]
```

 `ChipSum`在`examples/chipsumAI`中提供了基于全连接层的人工智能手写体识别示例，默认路径编译完成后，可以在`./build`中查看编译结果。
```
# default path
cd ./build
./ChipSumMnist
```
预期输出：
```
******input is****** : 9
 [                                                       ]
 [                                                       ]
 [                                                       ]
 [                                                       ]
 [                                                       ]
 [                                                       ]
 [                        # # # # # # # #                ]
 [                      # # # # # # # # # #              ]
 [                    # # # # # # # # # # #              ]
 [                    # # # # # # # # # # #              ]
 [                    # # # # #   # # # # #              ]
 [                    # # # # # # # # # # #              ]
 [                      # # # # # # # # # #              ]
 [                      # # # # # # # # # #              ]
 [                          # # # # # # #                ]
 [                          # # # # # #                  ]
 [                        # # # # # # #                  ]
 [                      # # # # # #                      ]
 [                    # # # # # #                        ]
 [                  # # # # # #                          ]
 [                # # # # # # #                          ]
 [                # # # # # #                            ]
 [              # # # # # #                              ]
 [            # # # # # #                                ]
 [              # # # # #                                ]
 [              # # # #                                  ]
 [                                                       ]
 [                                                       ]
*****prediction is***** : 9
```