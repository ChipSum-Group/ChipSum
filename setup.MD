setup.py编译方法：


例（成都超算）：
                
python3 setup.py hip=/public/software/compiler/rocm/rocm-3.9.1/  arch=VEGA906 compiler=/public/software/compiler/rocm/rocm-3.9.1/bin/hipcc
                
                
例（Volta72）:

python3 setup.py cuda=/Path/To/Your/Cuda arch=Volta72
                

例（Pascal61）@gassive：

python3 setup.py cuda=/Path/To/Your/Cuda arch=Pascal61
                
                
                
例（指定安装目录）：

mkdir /Some/Where/You/Like
python3 setup.py cuda=/Path/To/Your/Cuda arch=Pascal61 prefix=/Some/Where/You/Like



例（指定核编译，默认为4核）：
python3 setup.py cuda=/Path/To/Your/Cuda arch=Pascal61 j=32
