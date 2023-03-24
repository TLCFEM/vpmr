# VPMR C++ Implementation

## What Is This

This is a C++ implementation of the VPMR algorithm to compute the approximation of arbitrary smooth kernel.

Check the reference paper [10.1007/s10915-022-01999-1](https://doi.org/10.1007/s10915-022-01999-1) and
the [original](https://github.com/ZXGao97/VPMR) MATLAB implementation for more details.

## Dependency

1. [gmp](https://gmplib.org/) for multiple precision arithmetic
2. [mpfr](https://www.mpfr.org/) for multiple-precision floating-point computations
3. [mpreal](http://www.holoborodko.com/pavel/mpfr/) `mpreal` type C++ wrapper, included
4. [BigInt](https://github.com/faheel/BigInt) `BigInt` arbitrary large integer for combinatorial number, included
5. [Eigen](https://eigen.tuxfamily.org/) for matrix decomposition, included
6. [tbb](https://github.com/oneapi-src/oneTBB) for parallel computing

## How To Compile

The following is based on Fedora. For Windows users, please use WSL.
Native support is hard as porting dependencies to Windows is cumbersome.

1. Install compiler, `cmake` and `git`.
   ```bash
   sudo dnf install gcc g++ gfortran cmake git -y
   ```
2. Install `tbb`, `gmp` and `mpfr` packages.
   ```bash
   sudo dnf install tbb-devel mpfr-devel gmp-devel -y
   ```
3. Initialise submodules.
   ```bash
   git submodule update --init --recursive
   ```
4. Configure and compile.
   ```bash
   cmake -DCMAKE_BUILD_TYPE=Release .
   make
   ```

## How To Use

### Provide Kernel

It is necessary to compile the application with the desired kernel function first.

In the `Kernel` class (`src/VPMR.cpp`), change the `smooth_function` function to the target kernel function.

```cpp
class Kernel {
    //
    // ... other functions
    //
    
    // change this function to the target kernel function in convolution
    static mpreal smooth_function(const mpreal &x) {
        return exp(-x * x / 4);
    }
};
```

### Usage

All available options are:

```text
Usage: vpmr [options]

Options:

   -nc <int>    number of exponent (default: 4)
   -n <int>     number of terms (default: 10)
   -d <int>     number of digits (default: 512)
   -q <int>     quadrature order (default: 500)
   -e <float>   tolerance (default: 1E-8)
   -h           print this help message
```

### Example

The default kernel is `exp(-x^2/4)`. One can run the application with the following command:

```bash
./vpmr -n 30
```

The output is:

```text
Using the following parameters:
        nc = 4.
         n = 30.
     order = 500.
 precision = 512.
 tolerance = 1.0000e-08.

M = 
+1.0589202279681926e-11+0.0000000000000000e+00j
-5.4905134221689723e-03+2.2104939243740062e-03j
-5.4905134221689723e-03-2.2104939243740062e-03j
+1.1745193571738945e+01+2.3424418474362436e-153j
-5.5143304351134406e+00-5.7204056791636839e+00j
-5.5143304351134406e+00+5.7204056791636839e+00j
-1.6161617424833765e-02+2.3459542440459513e+00j
-1.6161617424833765e-02-2.3459542440459513e+00j
+1.6338578576177490e-01-1.9308431539218421e-01j
+1.6338578576177490e-01+1.9308431539218421e-01j
S = 
+0.0000000000000000e+00+0.0000000000000000e+00j
+1.7655956664692953e+00-2.7555720406099038e+00j
+1.7655956664692953e+00+2.7555720406099038e+00j
+1.8757961592204051e+00-0.0000000000000000e+00j
+1.8700580506914817e+00-6.2013413918954552e-01j
+1.8700580506914817e+00+6.2013413918954552e-01j
+1.8521958553280000e+00-1.2601975249082220e+00j
+1.8521958553280000e+00+1.2601975249082220e+00j
+1.8197653300065935e+00-1.9494562062795735e+00j
+1.8197653300065935e+00+1.9494562062795735e+00j

Running time: 11 s.
```

### Visualisation

The `plotter` folder contains a Python script to plot the result.

It is necessary to have `matplotlib` and `numpy` available in, for example, the virtual environment.

Copy and paste the `M` and `S` values to `weights` variable in `plotter/plotter.py` and run the script to plot the
result.

For the above example, the result is:

![exp(-x^2/4)](example.png)

## Binary

The binary released requires available `gmp`, `mpfr` and `tbb` libraries.

```bash
vpmr/cmake-build-release on  master [!] via △ v3.26.0 
❯ ldd vpmr 
        linux-vdso.so.1 (0x00007ffcf3121000)
        libgmp.so.10 => /lib64/libgmp.so.10 (0x00007f72087e8000)
        libmpfr.so.6 => /lib64/libmpfr.so.6 (0x00007f7208736000)
        libtbb.so.2 => /lib64/libtbb.so.2 (0x00007f72086f2000)
        libstdc++.so.6 => /lib64/libstdc++.so.6 (0x00007f7208400000)
        libm.so.6 => /lib64/libm.so.6 (0x00007f7208320000)
        libgcc_s.so.1 => /lib64/libgcc_s.so.1 (0x00007f72086d0000)
        libc.so.6 => /lib64/libc.so.6 (0x00007f7208143000)
        /lib64/ld-linux-x86-64.so.2 (0x00007f72088a1000)
```

You may wish to compile the application anyway since I do not know what kernel you want to approximate.
