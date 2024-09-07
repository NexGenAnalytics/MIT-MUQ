
import os
import sys
import shutil
import urllib.request
import subprocess
from typing import Union
from os import environ
from pathlib import Path
import concurrent.futures
import multiprocessing

def remove_everything_if_needed_from(pathdir, force_rebuild: bool):
    if os.path.exists(pathdir) and force_rebuild:
        print(f'removing all content in {pathdir} becuase -f was passed')
        shutil.rmtree(pathdir)

def get_full_path_to_cmake_config_dir(
    install_path_to_search: Union[str, os.PathLike],
    config_file_name: str
):
    for path in Path(install_path_to_search).rglob(config_file_name):
        return os.path.dirname(path)


def build_install_boost(
    workdir: str,
    force_rebuild: bool
):
    # boost is special since we cannot use cmake for this version

    myname = "boost"
    parent_path = os.path.join(workdir, myname)

    zip_name = "boost_1_85_0.tar.gz"
    url = f'https://archives.boost.io/release/1.85.0/source/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  "boost_1_85_0")
    install_path = os.path.join(parent_path, "install")

    # prepare
    remove_everything_if_needed_from(parent_path, force_rebuild)

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
        "--with-libraries=filesystem,system,graph")
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

    target = 'BoostConfig.cmake'
    return get_full_path_to_cmake_config_dir(install_path, target)


def build_install_impl(
        tplname: str,
        force_rebuild: bool,
        parentdir: Union[str, os.PathLike],
        url: str,
        zipname: Union[str, os.PathLike],
        zippath: Union[str, os.PathLike],
        unpacked: Union[str, os.PathLike],
        builddir: Union[str, os.PathLike],
        installdir: Union[str, os.PathLike],
        cmake_extra_args
):
    # prepare
    remove_everything_if_needed_from(parentdir, force_rebuild)

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
        "-DCMAKE_C_COMPILER=", os.environ['CC'],
        "-DCMAKE_CXX_COMPILER=", os.environ['CXX'],
        "-S", unpacked,
        "-B", builddir,
        f'-DCMAKE_INSTALL_PREFIX={installdir}')
    exeargs += cmake_extra_args
    print(exeargs)

    logfile_path = parentdir + "/logfile_build"
    logfile = open(logfile_path, "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    if p.returncode != 0:
        with open(logfile_path, 'r') as f:
            print(f.read())
    assert p.returncode == 0

    # make and install
    print("3. make and install")
    os.chdir(builddir)
    exeargs = ("make", "-j4", "install")
    logfile_path = parentdir + "/logfile_makeinstall"
    logfile = open(logfile_path, "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    if p.returncode != 0:
        with open(logfile_path, 'r') as f:
            print(f.read())
    assert p.returncode == 0

    print(f'success installing: {tplname}\n')


def build_install_hdf5(
    workdir: str,
    force_rebuild: bool
):
    myname = "hdf5"
    parent_path = os.path.join(workdir, myname)

    zip_name = "hdf5-1_14_1-2.zip"
    url = f'https://github.com/HDFGroup/hdf5/releases/download/hdf5-1_14_1-2/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  "hdfsrc")
    build_path = os.path.join(parent_path, "build")
    install_path = os.path.join(parent_path, "install")
    custom_cmake_args = (f'-DHDF5_BUILD_CPP_LIB=ON',)

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

    target = 'hdf5-config.cmake'
    return get_full_path_to_cmake_config_dir(install_path, target)


def build_install_nlopt(
    workdir: str,
    force_rebuild: bool
):
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

    target = 'NLoptConfig.cmake'
    return get_full_path_to_cmake_config_dir(install_path, target)


def build_install_sundials(
    workdir: str,
    force_rebuild: bool
):
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

    target = 'SUNDIALSConfig.cmake'
    return get_full_path_to_cmake_config_dir(install_path, target)


def build_install_eigen(
    workdir: str,
    force_rebuild: bool
):
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

    target = 'Eigen3Config.cmake'
    return get_full_path_to_cmake_config_dir(install_path, target)


def build_install_nanoflann(
    workdir: str,
    force_rebuild: bool
):
    myname = "nanoflann"
    parent_path = os.path.join(workdir, myname)

    zip_name = "v1.5.5.zip"
    url = f'https://github.com/jlblancoc/nanoflann/archive/refs/tags/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  "nanoflann-1.5.5")
    build_path = os.path.join(parent_path, "build")
    install_path = os.path.join(parent_path, "install")
    custom_cmake_args = (
        "-DNANOFLANN_BUILD_EXAMPLES=OFF",
        "-DNANOFLANN_BUILD_TESTS=OFF")

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

    target = 'nanoflannConfig.cmake'
    return get_full_path_to_cmake_config_dir(install_path, target)


def build_install_parcer(
    workdir: str,
    force_rebuild: bool
):
    myname = "parcer"
    parent_path = os.path.join(workdir, myname)

    commit = "3b1ee6dc3d73"
    zip_name = f'{commit}.zip'
    url = f'https://bitbucket.org/mituq/parcer/get/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  f'mituq-parcer-{commit}')
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

    target = 'PARCERConfig.cmake'
    return get_full_path_to_cmake_config_dir(install_path, target)



def build_install_spdlog(
    workdir: str,
    force_rebuild: bool
):
    myname = "spdlog"
    parent_path = os.path.join(workdir, myname)

    zip_name = f'v1.10.0.zip'
    url = f'https://github.com/gabime/spdlog/archive/refs/tags/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  f'spdlog-1.10.0')
    build_path = os.path.join(parent_path, "build")
    install_path = os.path.join(parent_path, "install")
    custom_cmake_args = (
        '-DCMAKE_POSITION_INDEPENDENT_CODE=ON',
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

    target = 'spdlogConfig.cmake'
    return get_full_path_to_cmake_config_dir(install_path, target)




def build_install_otf2(
    workdir: str,
    force_rebuild: bool
):
    myname = "otf2"
    parent_path = os.path.join(workdir, myname)

    zip_name = f'otf2-3.0.3.tar.gz'
    url = f'https://perftools.pages.jsc.fz-juelich.de/cicd/otf2/tags/otf2-3.0.3/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  f'otf2-3.0.3')
    build_path = os.path.join(parent_path, "build")
    install_path = os.path.join(parent_path, "install")

    # prepare
    remove_everything_if_needed_from(parent_path, force_rebuild)

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
    # note that we disable Python support here
    os.environ['PYTHON'] = ":"
    exeargs = ("./configure", 
        f'--prefix={install_path}', 
        f'CC={os.environ['CC']}', 
        f'CXX={os.environ['CXX']}')
    print(exeargs)

    logfile = open(parent_path + "/logfile_config", "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    assert p.returncode == 0

    # make and install
    print("3. make and install")
    exeargs = ("make", "-j4", "install")
    logfile = open(parent_path + "/logfile_makeinstall", "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    assert p.returncode == 0

    print(f'success installing: {myname}\n')

    return install_path



def build_install_stanmath(
        workdir: str,
        force_rebuild: bool
):
    myname = "stanmath"
    parent_path = os.path.join(workdir, myname)

    zip_name = "v2.18.0.zip"
    url = f'https://github.com/stan-dev/math/archive/refs/tags/{zip_name}'
    zip_path = os.path.join(parent_path, zip_name)
    unpack_path = os.path.join(parent_path,  "math-2.18.0")

    # prepare
    remove_everything_if_needed_from(parent_path, force_rebuild)

    if os.path.exists(parent_path):
        print(f'skipping {parent_path}, because found already. Use -f to rebuid.')

    print("-"*50)
    print(f'WORKING ON: {myname}')
    print("-"*50)

    os.mkdir(parent_path)

    # fetch
    print("1. fetching")
    urllib.request.urlretrieve(url, zip_path)

    # unpacking
    print("2. unpacking")
    shutil.unpack_archive(zip_path, parent_path)

    print(f'success getting: {myname}\n')

    print(unpack_path)
    return os.path.abspath(unpack_path)


def check_compilers(compiler_type: str):
    compiler_alias_found = environ.get(compiler_type)
    if compiler_alias_found is not None:
        print(f'found {compiler_alias_found}')
        return

    if not compiler_alias_found:
        print(f'-'*50)
        print(f'*** FATAL ERROR ***')
        print(f'-'*50)
        print(f'{compiler_type} not found in environment, I cannot proceed!')
        if compiler_type == 'CC':
            print(f'please set {compiler_type} to a valid C compiler')
        elif compiler_type == 'CXX':
            print(f'please set {compiler_type} to a valid C++ compiler')
        else:
            print(f'Specified compiler type "{compiler_type}" is invalid. Use "CC" or "CXX".')
        print(f'-'*50)
        exit(11)


currently_supported_tpls = [
    "hdf5", "nlopt", "boost", "sundials", "eigen",
    "nanoflann", "stanmath", "parcer", "spdlog", "otf2"
]

def add_arguments(parser):
    parser.add_argument(
        "--wdir", "--workdir",
        dest="workdir",
        required=True,
        help="Full path to your desired working directory."
        "If the directory does not exist, it is created.",
    )

    parser.add_argument(
        "-f",
        dest="force_rebuild",
        required=False,
        action="store_true",
        help="Use it to delete everything inside the working directory, "
        "and build everything from scratch."
    )

    parser.add_argument(
        "--with",
        dest="target_list",
        type=str,
        nargs="*",
        required=True,
        default=["all"],
        choices=currently_supported_tpls + ['all'],
        help="List of libraries to build/install. Default = all. \n"
        "Use '--with all' to build everything.",
    )


def dispatch(
    workdir: Union[str, os.PathLike], 
    force_rebuild: bool, 
    tpl: str
):
    if tpl == 'hdf5':
        p = build_install_hdf5(workdir, force_rebuild)
        return ['HDF5_DIR', p, tpl + ' path']

    if tpl == 'nlopt':
        p = build_install_nlopt(workdir, force_rebuild)
        return ['NLopt_DIR', p, tpl + ' path']

    if tpl == 'boost':
        p = build_install_boost(workdir, force_rebuild)
        return ['Boost_DIR', p, tpl + ' path']

    if tpl == 'sundials':
        p = build_install_sundials(workdir, force_rebuild)
        return ['SUNDIALS_DIR', p, tpl + ' path']

    if tpl == 'eigen':
        p = build_install_eigen(workdir, force_rebuild)
        return ['Eigen3_DIR', p, tpl + ' path']

    if tpl == 'nanoflann':
        p = build_install_nanoflann(workdir, force_rebuild)
        return ['nanoflann_DIR', p, tpl + ' path']

    if tpl == 'stanmath':
        p = build_install_stanmath(workdir, force_rebuild)
        return ['stanmath_SRC_DIR', p, tpl + ' path']

    if tpl == 'parcer':
        p = build_install_parcer(workdir, force_rebuild)
        return ['PARCER_DIR', p, tpl + ' path']

    if tpl == 'spdlog':
        p = build_install_spdlog(workdir, force_rebuild)
        return ['spdlog_DIR', p, tpl + ' path']

    if tpl == 'otf2':
        p = build_install_otf2(workdir, force_rebuild)
        return ['otf2_DIR', p, tpl + ' path']


##########################################################
if __name__ == "__main__":
##########################################################
    import argparse

    parser = argparse.ArgumentParser(
        prog="python bla bla"
    )
    add_arguments(parser)
    args = parser.parse_args()

    #
    # check things
    assert os.path.isabs(args.workdir)
    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)
    check_compilers('CC')
    check_compilers('CXX')

    #
    # summary
    final_tpls = currently_supported_tpls if 'all' in args.target_list \
        else args.target_list
    print(f'workdir set to {args.workdir}')
    print(f'TPLs to build {final_tpls}')

    pool_size = 4
    with concurrent.futures.ProcessPoolExecutor(max_workers = pool_size) as executor:
            these_futures = [executor.submit(dispatch, args.workdir, args.force_rebuild, tpl) for tpl in final_tpls]
            concurrent.futures.wait(these_futures)

    results = [future.result() for future in these_futures]
    print(results)


    paths_for_config = {}
    for it in results:
        paths_for_config[it[0]] = [it[2], it[1]]

    print(f'-'*50)
    print(f'writing cmake cache file for selected tpls')
    cache_file_out = os.path.join(args.workdir, 'tpls_cache.txt')
    with open(cache_file_out, 'w') as f:
        for k,v in paths_for_config.items():
            f.write(f'set({k} {v[1]} CACHE PATH "{v[0]}")\n')
    print(f'-'*50)
    print(f'use "-C {cache_file_out}" to find the selected TPLs when building MUQ')
