
import os
import shutil
import urllib.request
import subprocess


def rmeverything_if_needed(pathdir, force_rebuild: bool):
    if os.path.exists(pathdir) and force_rebuild:
        print(f'removing all content in {pathdir} becuase -f was passed')
        shutil.rmtree(pathdir)
        
def build_install_impl(
        tplname: str,
        force_rebuild: bool,
        parentdir: str | os.PathLike,
        url: str,
        zipname: str | os.PathLike,
        zippath: str | os.PathLike,
        unpacked: str | os.PathLike,
        builddir: str | os.PathLike,
        installdir: str | os.PathLike,
        cmake_extra_args
):
    # prepare
    rmeverything_if_needed(parentdir, force_rebuild)

    if os.path.exists(parentdir):
        print(f'skipping {parentdir}, because found already. Use -f to rebuid.')

    print("-"*50)
    print(f'WORKING ON: {tplname}')
    print("-"*50)

    os.mkdir(parentdir)

    # fetch
    print("1. fetching")
    urllib.request.urlretrieve(url, zippath)
    shutil.unpack_archive(zippath, parentdir)

    # configure
    print("2. configuring")        
    exeargs = (
        "cmake",
        "-S", unpacked, 
        "-B", builddir,
        f'-DCMAKE_INSTALL_PREFIX={installdir}')
    exeargs += cmake_extra_args
    print(exeargs)
    
    logfile = open(parentdir + "/logfile_build", "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    assert p.returncode == 0

    # make and install
    print("3. make and install")            
    os.chdir(builddir)
    exeargs = ("make", "-j4", "install")
    logfile = open(parentdir + "/logfile_makeinstall", "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    assert p.returncode == 0

    print(f'success installing: {tplname}\n')
    

def build_install_hdf5(workdir: str, force_rebuild: bool):
    myname = "hdf5"    
    parent_path = os.path.join(workdir, myname)
    
    zip_name = "hdf5-1_8_19.zip"
    url = f'https://github.com/HDFGroup/hdf5/archive/refs/tags/{zip_name}'        
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  "hdf5-hdf5-1_8_19")    
    build_path = os.path.join(parent_path, "build")
    install_path = os.path.join(parent_path, "install")
    custom_cmake_args = (f'-DHDF5_BUILD_CPP_LIB=ON',)
    
    build_install_impl(
        force_rebuild,
        parent_path,
        url,
        zip_name,
        zip_path,
        unpack_path,
        build_path,
        install_path,
        custom_cmake_args
    )
    
def build_install_nlopt(workdir: str, force_rebuild: bool):
    myname = "nlopt"
    parent_path = os.path.join(workdir, myname)
    
    zip_name = "v2.8.0.zip"
    url = f'https://github.com/stevengj/nlopt/archive/refs/tags/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  "nlopt-2.8.0")    
    build_path = os.path.join(parent_path, "build")
    install_path = os.path.join(parent_path, "install")
    custom_cmake_args = ()
    
    build_install_impl(
        myname,
        force_rebuild,
        parent_path,
        url,
        zip_name,
        zip_path,
        unpack_path,
        build_path,
        install_path,
        custom_cmake_args
    )

def build_install_sundials(workdir: str, force_rebuild: bool):
    myname = "sundials"
    parent_path = os.path.join(workdir, myname)
    
    zip_name = "v5.5.0.zip"
    url = f'https://github.com/LLNL/sundials/archive/refs/tags/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  "sundials-5.5.0")    
    build_path = os.path.join(parent_path, "build")
    install_path = os.path.join(parent_path, "install")
    custom_cmake_args = (
        "-DBUILD_CVODE=OFF", 
        "-DBUILD_CVODES=ON",
        "-DBUILD_IDA=OFF",
        "-DBUILD_IDAS=ON",
        "-DBUILD_KINSOL=ON",
    )
    
    build_install_impl(
        myname,
        force_rebuild,
        parent_path,
        url,
        zip_name,
        zip_path,
        unpack_path,
        build_path,
        install_path,
        custom_cmake_args
    )

def build_install_eigen(workdir: str, force_rebuild: bool):
    myname = "eigen"
    parent_path = os.path.join(workdir, myname)
    
    zip_name = "eigen-3.3.7.zip"
    url = f'https://gitlab.com/libeigen/eigen/-/archive/3.3.7/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  "eigen-3.3.7")
    build_path = os.path.join(parent_path, "build")
    install_path = os.path.join(parent_path, "install")
    custom_cmake_args = ()
    
    build_install_impl(
        myname,
        force_rebuild,
        parent_path,
        url,
        zip_name,
        zip_path,
        unpack_path,
        build_path,
        install_path,
        custom_cmake_args
    )
    
    
def build_install_boost(workdir: str, force_rebuild: bool):
    # boost is special since we cannot use cmake for this version
    
    myname = "boost"
    parent_path = os.path.join(workdir, myname)
    
    zip_name = "boost_1_85_0.tar.gz"
    url = f'https://archives.boost.io/release/1.85.0/source/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  "boost_1_85_0")
    install_path = os.path.join(parent_path, "install")

    # prepare
    rmeverything_if_needed(parent_path, force_rebuild)

    if os.path.exists(parent_path):
        print(f'skipping {parent_path}, because found already. Use -f to rebuid.')

    print("-"*50)
    print(f'WORKING ON: {myname}')
    print("-"*50)

    os.mkdir(parent_path)

    # fetch
    print("1. fetching")
    urllib.request.urlretrieve(url, zip_path)
    shutil.unpack_archive(zip_path, parent_path)

    os.chdir(unpack_path)
    
    # configure
    print("2. configuring")
    exeargs = (
        "./bootstrap.sh",
        f'--prefix={install_path}',
        "--with-libraries=graph")
    print(exeargs)
    
    logfile = open(parent_path + "/logfile_config", "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    assert p.returncode == 0

    # make and install
    print("3. make and install")            
    exeargs = ("./b2", "-j4", "install")
    logfile = open(parent_path + "/logfile_makeinstall", "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    assert p.returncode == 0

    print(f'success installing: {myname}\n')
    

def add_arguments(parser):
    parser.add_argument(
        "--wdir",
        dest="workdir",
        required=True,
        help="Full path to your desired working directory."
        "If the directory does not exist, it is created.",
    )
    parser.add_argument(
        "-f", dest="force_rebuild", required=False, action="store_true"
    )
    
    parser.add_argument(
        "--with-hdf5", dest="with_hdf5", required=False, action="store_true"
    )
    
    parser.add_argument(
        "--with-nlopt", dest="with_nlopt", required=False,action="store_true"
    )
    
    parser.add_argument(
        "--with-boost", dest="with_boost", required=False,action="store_true"
    )
    
    parser.add_argument(
        "--with-sundials", dest="with_sundials", required=False,action="store_true"
    )

    parser.add_argument(
        "--with-eigen", dest="with_eigen", required=False,action="store_true"
    )
    

##########################################################
if __name__ == "__main__":
##########################################################
    import argparse

    parser = argparse.ArgumentParser(
        prog="python bla bla"
    )
    add_arguments(parser)
    args = parser.parse_args()

    # check things
    assert os.path.isabs(args.workdir)
    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)

    print(args.workdir)
    if (args.with_hdf5):
        build_install_hdf5(args.workdir, args.force_rebuild)

    if (args.with_nlopt):
        build_install_nlopt(args.workdir, args.force_rebuild)        

    if (args.with_boost):
        build_install_boost(args.workdir, args.force_rebuild)

    if (args.with_sundials):
        build_install_sundials(args.workdir, args.force_rebuild)                

    if (args.with_eigen):
        build_install_eigen(args.workdir, args.force_rebuild)                
