
import os
import sys


def print_help():
    with open("./setup.MD") as f:
        print(f.read())
    
    exit(0)



if sys.argv[1] in ["--help","-h"]:
    print_help()





cuda_path = None
hip_path = None
arch = None
chipsum_prefix = None
chipsum_compiler = None
make_procs = 4





for i in range(1,len(sys.argv)):
    key,val=sys.argv[i].split("=")

    if key.lower() == "cuda":
        cuda_path = val
    elif key.lower() == "prefix":
        chipsum_prefix = val
    elif key.lower()=="hip":
        hip_path = val
    elif key.lower()=="compiler":
        chipsum_compiler = val
    elif key.lower()=="j":
    	make_procs = int(val)
    elif key.lower()=="prefix":
    	chipsum_prefix=val
    if key.lower()=="arch":
        arch=val
    






def kokkos_source_build():
    if os.path.exists("./tpl/kokkos/build"):
    	os.system("rm -rf ./tpl/kokkos/build")
    	
    os.mkdir("./tpl/kokkos/build")
    if os.path.exists("./tpl/kokkos-build"):
    	os.system("rm -rf ./tpl/kokkos-build")    
    os.mkdir("./tpl/kokkos-build")
        
 
    org_path = os.path.abspath(".")
    os.chdir("./tpl/kokkos/build")
    
    bash = "../generate_makefile.bash "
      
    root = "--kokkos-path=../ "
    
    prefix = "--prefix=../../kokkos-build "
    
    with_gpu = ""
    
    if cuda_path != None:
        with_gpu = "--with-cuda="+cuda_path+" "
    elif hip_path != None:
        with_gpu = "--with-hip="+hip_path+" "
    
    kokkos_arch = "--arch="+arch+" "
    
    if cuda_path:
        compiler = "--compiler="+org_path+"/tpl/kokkos/bin/nvcc_wrapper "
    else:
        compiler = "--compiler="+chipsum_compiler+" "
    
    os.system(bash+root+prefix+with_gpu+kokkos_arch+compiler+"--cxxstandard=14 --with-openmp --with-serial --disable-tests")
    os.system("make -j"+str(make_procs))
    os.system("make install")
    os.chdir(org_path)
    
    
    
def kokkoskernels_source_build():
    if os.path.exists("./tpl/kokkos-kernels/build"):
    	os.system("rm -rf ./tpl/kokkos-kernels/build")
    os.mkdir("./tpl/kokkos-kernels/build")

    if os.path.exists("./tpl/kokkos-kernels-build"):
    	os.system("rm -rf ./tpl/kokkos-kernels-build")
    os.mkdir("./tpl/kokkos-kernels-build")    
    org_path = os.path.abspath(".")
    os.chdir("./tpl/kokkos-kernels/build")
    
    bash = "../cm_generate_makefile.bash "

    
    
    kokkos_prefix = "../../kokkos-build "
    
    prefix = "--prefix=../../kokkos-kernels-build "
    
    with_gpu = ""
    
    if cuda_path != None:
        with_gpu = "--with-cuda="+cuda_path+" "
    elif hip_path != None:
        with_gpu = "--with-hip="+hip_path+" "
   
    if arch!=None:
        kokkos_arch = "--arch="+arch+" "
    else:
        kokkos_arch = ""

    compiler = ""        
    if cuda_path:
        compiler = "--compiler="+org_path+"/tpl/kokkos/bin/nvcc_wrapper "
    else:
        compiler = "--compiler="+chipsum_compiler+" "
    
    os.system(bash+kokkos_prefix+prefix+with_gpu+kokkos_arch+compiler+"--cxxstandard=14 --with-openmp --with-serial --release --kokkoskernels-path=../ --disable-tests")
    os.system("make -j"+str(make_procs))
    os.system("make install")
    os.chdir(org_path)
    
tpl_flags = [
        "kokkos",
        "kokkoskernels"
        ]
    
tpl_build_flags = [
        kokkos_source_build,
        kokkoskernels_source_build
        ]

tpl_source_flags = [
        "./tpl/kokkos/CMakeLists.txt",
        "./tpl/kokkos-kernels/CMakeLists.txt"
        ]


tpl_lib_flags = [
        "./tpl/kokkos-build/lib/cmake/Kokkos/KokkosConfig.cmake",
        "./tpl/kokkos-kernels-build/lib/cmake/KokkosKernels/KokkosKernelsConfig.cmake"
        ]
        
tpl_lib64_flags = [
        "./tpl/kokkos-build/lib64/cmake/Kokkos/KokkosConfig.cmake",
        "./tpl/kokkos-kernels-build/lib64/cmake/KokkosKernels/KokkosKernelsConfig.cmake"
        ]



    
    
    
def build_all():
    
    for i in range(len(tpl_lib_flags)):
        if os.path.exists(tpl_lib_flags[i]) or os.path.exists(tpl_lib64_flags[i]):
            pass
        else:
            
            if os.path.exists(tpl_source_flags[i]):
                tpl_build_flags[i]()
            
            else:
                print("You don't have "+tpl_flags[i]+" in ChipSum/tpl directory.")
                exit(-1)
        

        
    
    if os.path.exists("build"):
       os.system("rm -rf build")
    os.mkdir("build")
    org_path = os.path.abspath(".")
    
    build_prefix = ""
    
    if chipsum_prefix!=None:
        build_prefix = "-DCMAKE_INSTALL_PREFIX="+chipsum_prefix+" "
   
    if os.path.exists(tpl_lib_flags[0]): 
        os.chdir("build")
        os.system("cmake -DChipSum_ENABLE_KokkosKernels=yes "+ build_prefix+"..")
    elif os.path.exists(tpl_lib64_flags[0]): 
        os.chdir("build")
        os.system("cmake -DChipSum_ENABLE_KokkosKernels64=yes "+ build_prefix+"..")
    os.system("make -j"+str(make_procs))
    if build_prefix != "":
        os.system("make install")
    os.chdir(org_path)
        
        
        
            
build_all()

