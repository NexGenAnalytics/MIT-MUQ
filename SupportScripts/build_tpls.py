
import os
import shutil
import urllib.request
import subprocess

def add_arguments(parser):
    parser.add_argument(
        "--wdir",
        dest="workdir",
        required=True,
        help="Full path to your desired working directory."
        "If the directory does not exist, it is created.",
    )

    parser.add_argument(
        "--with-hdf5",
        dest="with_hdf5",
        required=False,
        action="store_true",
        help="Whether to build and install HDF5.",
    )

    parser.add_argument(
        "-f",
        dest="force_rebuild",
        required=False,
        action="store_true"
    )

def rmeverything_if_needed(pathdir, force_rebuild: bool):
    if os.path.exists(pathdir) and force_rebuild:
        print(f'removing all content in {pathdir} becuase -f was passed')
        shutil.rmtree(pathdir)
        
def build_install_hdf5(workdir: str, force_rebuild: bool):
    hdf5_parentdir = os.path.join(workdir, "hdf5")
    zipname = "hdf5-1_8_19.zip"
    hdf5_zippath = hdf5_parentdir + "/" + zipname
    hdf5_unpacked = os.path.join(hdf5_parentdir,  "hdf5-hdf5-1_8_19")    
    hdf5_builddir = os.path.join(hdf5_parentdir, "build")
    hdf5_installdir = os.path.join(hdf5_parentdir, "install")
    
    # prepare
    rmeverything_if_needed(hdf5_parentdir, force_rebuild)

    if os.path.exists(hdf5_parentdir):
        print(f'skipping hdf5, because found already. Use -f to rebuid.')
        
    os.mkdir(hdf5_parentdir)

    # fetch
    url = f'https://github.com/HDFGroup/hdf5/archive/refs/tags/{zipname}'
    urllib.request.urlretrieve(url, hdf5_zippath)
    shutil.unpack_archive(hdf5_zippath, hdf5_parentdir)

    # configure 
    exeargs = (
        "cmake",
        "-S", hdf5_unpacked, 
        "-B", hdf5_builddir,
        f'-DCMAKE_INSTALL_PREFIX={hdf5_installdir}',
        f'-DHDF5_BUILD_CPP_LIB=ON')
    print(exeargs)
    logfile = open(hdf5_parentdir + "/logfile_build", "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    assert p.returncode == 0

    # make and install
    os.chdir(hdf5_builddir)
    exeargs = ("make", "-j4", "install")
    logfile = open(hdf5_parentdir + "/logfile_makeinstall", "w")
    p = subprocess.Popen(exeargs, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
    assert p.returncode == 0

    print("hdf5 success")
    
    
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
    
