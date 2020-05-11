# Contents
* [Cloning NL-PartSol](#cloning-NL-PartSol)
* [NL-PartSol Dependencies](#NL-PartSol-dependencies)
  * [Linux Installation](#linux-installation)
  * [Windows Installation](#windows-installation)
* [Basic Configuration](#basic-configuration)
* [Adding Applications](#adding-applications)
* [Adding NL-PartSol to Path](#adding-NL-PartSol-to-path)
* [Examples](#examples)
  * [Linux](#linux)
  * [Windows](#windows)
  * [MacOS](#macos)

## Cloning NL-PartSol

In order to obtain the source code of NL-PartSol, you will need to clone the repository using git.

		You can install git through the following command in Linux:
```Shell
sudo apt-get install git
```
In Windows, you can download it in:

* [Download Git](https://git-scm.com/downloads)


		Once git is installed you can fetch the code by using this command in a terminal:

		```Shell
git clone https://github.com/migmolper/NL-PartSol.git
```

## NL-PartSol Dependencies
	These are the basic dependecies needed to compile the Kratos Core and most of the Kratos applications.
		* C compiler
		* CMake
	 
	Additionaly, Visual Studio is required to compile in Windows.
		
- #### Linux installation

	The command below will install all the packages needed.

	```Shell
		sudo apt-get install gcc cmake
	```

	- #### Windows installation

	- Visual Studio

	*Visual Studio* is the only compiler officially supported to build *Kratos* under *Windows*. 
	The minimium required version is Visual Studio 2017, but we recommend to use Visual Studio 2019 or higher. 

	* [Download Visual Studio](https://visualstudio.microsoft.com/en/thank-you-downloading-visual-studio/?sku=Community&rel=15)

	Since *Visual Studio* is a multi-language IDE, some distributions come without C++ compiler. 
	Please, make sure that you can create a C++ project before continuing, in case C++ packages were missing you will be prompt to download them.

	- CMake
   * [Download CMake](http://cmake.org/download/)

		Once installing, please <span style="color:red"> do not forget to mark the option: '''"Add CMake to the system PATH for all users"'''</span>
   
		Minimum required version: CMake 3.14

	
## Basic Configuration

		You can find the new kratos configuration file in Kratos `scripts` folder: `standard_configure.sh` for linux, `standard_configure_max.sh` for MacOS, `standard_configure.bat` for win and others.

		Out of the box Kratos will try to find all necessary libraries in your system automatically, but we recommend you to copy these scripts and modify it according to your preferences. 
		Please take a look at the following configuration options:

		`KRATOS_BUILD_TYPE`

		Compilation Type. Options are `Release`,`RelWithDebInfo`,`Debug`,`FullDebug`,`Custom`

		**Release**: Full Release with maximum optimization options and no debug Info.

		**RelWithDebInfo**: Full Release with optimization and debug info. Adecuate to debug simple problems without losing performance.

		**Debug**: Debug build with no optimization flags.

		**FullDebug**: Debug build with no optimization falgs, extended debug info and extremly low performance.

		**Custom**: No flags are automatically added.

		`PYTHON_EXECUTABLE`

		Path to the python executable that Kratos will use. We recommend that you manually set this in case you have multiple versions of python in the system.
		Ubuntu users need to be extra careful with this as default versions tends to be Python2, while Kratos is preferably compiled with Python3

		`BOOST_ROOT`

		Path to boost root directory.


## Adding NL-PartSol to Path

As Kratos is not an executable but a set of modules and libraries, you will need to add them to the path. 
In order to do that please add the Kratos install folder (If you didn't touch anything should be `$KRATOS_SOURCE/bin/Release`)

```bash
export PYTHONPATH=$PYTHONPATH:$HOME/Kratos/bin/Release
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Kratos/bin/Release/libs
```

If you are in windows instead do:

```cmd
set PYTHONPATH=%PYTHONPATH%;C:/Kratos/bin/Release
set PATH=%PATH%;C:/Kratos/bin/Release/libs
```

## Examples
	
These examples are also located [in the /examples folder](https://github.com/migmolper/NL-PartSol/tree/master/examples). You can simply create your own copy:

```Shell
cp /path_to_NL_PartSol/scripts/standard_configure.sh /path_to_kratos/scripts/configure.sh
```

### Linux

```bash
# Function to add apps
add_app () {
export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
	}

# Set compiler
export CC=gcc
export CXX=g++

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE="Release"
export PYTHON_EXECUTABLE="/usr/bin/python3"

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/EigenSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" -DUSE_MPI=OFF -DUSE_EIGEN_MKL=OFF

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j4

```

### Windows

```cmd
rem Set compiler
set CC=cl.exe
set CXX=cl.exe

rem Set variables
set KRATOS_SOURCE=~0,-1%/..
set KRATOS_BUILD=%KRATOS_SOURCE%/build
set KRATOS_APP_DIR=applications

rem Set basic configuration
set KRATOS_BUILD_TYPE=Release
set BOOST_ROOT=C:\boost_1_67_0
set PYTHON_EXECUTABLE=C:\Python37\python.exe

rem Set applications to compile
set KRATOS_APPLICATIONS=
CALL :add_app %KRATOS_APP_DIR%\EigenSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\FluidDynamicsApplication;

rem Clean
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

rem Configure
@echo on
cmake -G"Visual Studio 16 2019" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"  ^
-DUSE_EIGEN_MKL=OFF

rem Build
cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64
goto:eof

rem Function to add apps
:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof

```

### MacOS

```bash
# Function to add apps
add_app () {
export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
	}

# Set compiler
export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE="Release"
export BOOST_ROOT="/path/to/boost"
export PYTHON_EXECUTABLE="/Library/Frameworks/Python.framework/Versions/3.7/bin/python3"

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/EigenSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
/Applications/CMake.app/Contents/bin/cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 -L/usr/local/opt/llvm/lib" \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 -L/usr/local/opt/llvm/lib" \
-DUSE_EIGEN_MKL=OFF

# Buid
/Applications/CMake.app/Contents/bin/cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j3

```
