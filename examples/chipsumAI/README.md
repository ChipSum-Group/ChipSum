
##编译使用mnist

### 编译`ChipSum`

首先需要将`ChipSum`编译并make install至指定文件夹，如：examples/chipsumAI/install_lib

```
cd path_to_chipsum
cd ./build
cmake -DCMAKE_INSTALL_PREFIX=../examples/chipsumAI/install_lib ..
make install
ll ../examples/chipsumAI/install_lib
```
成功make install后，会在目标文件夹下(例子中，install_lib文件夹下)新生成以下文件。

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
此时，若需在其他文件中使用`ChipSum`，仅需`#include "ChipSum.hpp"`即可使用相应数据结构和函数。


### 编译`chipsumAI`
在运行手写体识别代码样例前，需要对项目进行编译，已获得可运行的可执行文件。

```
cd path_to_chipsum
cd ./examples/chipsumAI
mkdir build && cd build

# ChipSum_DIR should be absolute path
export ChipSum_DIR=/path/to/chipsum/

# ChipSumLib_DIR should be absolute path
export ChipSumLib_DIR=/path/to/install_lib/

cmake -DChipSum_DIR=${ChipSum_DIR} -DChipSumLib_DIR=${ChipSumLib_DIR} ..

make -j8

./mnist
```

使用命令`./mnist`运行mnist实现手写体识别，预期输出示例如下：

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