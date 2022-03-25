


##编译使用mnist

### 编译ChipSum

- `chipsumAI`依赖于`ChipSum`的数据结构等，故需要先编译`ChipSum`，并make install 到指定文件夹，如：ChipSum文件夹内。
```
# default path
cd ./build
cmake -DCMAKE_INSTALL_PREFIX=../ ..
make install
```
成功make install后，会在目标文件夹下新生成以下文件。
```
.
├── bin
│   └── ChipSum
├── include
│   ├── chipsum
│   └── ChipSum.hpp
├── lib
│   └── libchipsum.a
```
在其他文件中使`ChipSum`，仅需`#include "ChipSum.hpp"`即可使用相应数据结构和函数。


### 编译chipsumAI
- 在运行手写体识别代码样例前，需要对项目进行编译，已获得可运行的可执行文件。
```
# default path
cd ./examples/chipsumAI
mkdir build && cd build

export ChipSum_DIR=/path/you/want/

cmake -DKokkos_DIR=${ChipSum_DIR}/tpls/kokkos-build/lib/cmake/Kokkos -DKokkosKernels_DIR=${ChipSum_DIR}/tpls/kokkos-kernels-build/lib/cmake/KokkosKernels -DChipSum_DIR=${ChipSum_DIR} ..

make -j8
./mnist
```
运行mnist实现手写体识别，预期输出如：
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