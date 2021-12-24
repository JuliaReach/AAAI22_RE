# Sherlock experiments

This manual describes the steps to run the Sherlock experiments.



## Installation

- Install the [Gurobi solver](https://www.gurobi.com/) on your system.

- Install GSL. On Debian (Linux) systems the command is:

```bash
sudo apt-get install libgsl-dev
```

- Install GLPK. On Debian (Linux) systems the command is:

```bash
sudo apt-get install libglpk-dev
```

- Clone the repository from [here](https://bitbucket.org/souradeep-111/hscc_2019_implementation/src/master/).

```bash
git clone https://bitbucket.org/souradeep-111/hscc_2019_implementation.git
```

- Copy the folders from where the repository containing this manual to the new Sherlock folder and overwrite the files if necessary.

- Go to the folder `neural_network_reachability/examples/tests/`.

- Modify `Makefile.locale` containing information to find your Gurobi installation. We used Gurobi 8.0.1 because that version was also used by the Sherlock authors. Other versions may also be compatible. Example:

```bash
GUROBI_PATH = $HOME/programs/gurobi-8.0.1
HOST_ARCH = linux64
```

- The following steps are already done if you copied over the files in the step above.

    - Modify `Makefile` by changing the Gurobi paths and adding the `-no-pie` option (see also the _Troubleshooting_ section below):

    ```bash
    +include Makefile.locale
    +GUROBI_INCLUDEDIR=$(strip $(GUROBI_PATH))/$(strip $(HOST_ARCH))/include/
    +GUROBI_LIBDIR=$(strip $(GUROBI_PATH))/$(strip $(HOST_ARCH))/lib/

    LIBS = -lgurobi_c++ -lgurobi80 -lflowstar -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk -lmpfi -D_GLIBCXX_USE_CXX11_ABI=0 -m64 -w
    CFLAGS = -I . -I ./headers -I /usr/local/include/ -I $(GUROBI_INCLUDEDIR) \
        -I ./headers/eigen_file/ -I ./headers/neural_rule_analysis/src_new -g -O3 -std=c++11
    LINK_FLAGS = -g -L ./ -L /usr/local/lib/ -L $(GUROBI_LIBDIR) -L ../flowstar-release -no-pie
    ```
    - Download `Eigen` from [here](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) and extract the archive into the folder `headers/eigen_file`.

    - Add a file `eigen_ridge.hpp` with the following content to the folder `headers` (see [here](https://github.com/souradeep-111/sherlock/issues/5#issuecomment-649077490)):

    ```c
    #ifndef _EIGEN_RIDGE_HPP_

    #define _EIGEN_RIDGE_HPP_

    /*
    * L2-regularized (ridge) linear regression without intercept in c++11 as Eigen3 template function
    * Uses singular value decomposition, works with both tall and wide design matrices
    * Written by Carlos Guerreiro carlos@perceptiveconstructs.com
    * This is free and unencumbered software released into the public domain.
    */

    #include <Eigen/Dense>

    template<typename M, typename V, typename P>
    M ridge(const M& A, const V& y,  P alpha) {
    const auto& svd = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    const auto& s = svd.singularValues();
    const auto r = s.rows();
    const auto& D = s.cwiseQuotient((s.array().square() + alpha).matrix()).asDiagonal();
    return svd.matrixV().leftCols(r) * D * svd.matrixU().transpose().topRows(r) * y;
    }

    #endif
    ```
    - Fix the file `src/RangeToVariables.cpp`:

       - Add `#include <cassert>`.
       - Replace `string s;` by `std::string s;`.

- Run `make benchmarks`.



## Troubleshooting

- If the `make` process fails with errors of the form "undefined reference to `GRBModel::`", one needs to build Gurobi from source (see [here](https://github.com/souradeep-111/sherlock/issues/3), which points [here](https://groups.google.com/g/gurobi/c/9RXVpObMJxM)). For this, go to the Gurobi folder and there to `src/build/`, run `make`, and copy the resulting `libgurobi_c++.a` to Gurobi's `lib` folder (after renaming the old file).

- If the `make` process fails with errors of the form "can not be used when making a PIE object; recompile with -fPIE", add the option `-no-pie` to the linking command in `Makefile` under `LINK_FLAGS` (so it reads `LINK_FLAGS = -g -L ./ -L /usr/local/lib/ -L $(GUROBI_LIBDIR) -L ./flowstar-release/ -no-pie`).



## Running the experiments

Sherlock must be run with a pointer to Gurobi. One way to do that is to add Gurobi to `LD_LIBRARY_PATH`. Another way is to provide the path to the run script, e.g.:

```bash
LD_LIBRARY_PATH=link_to_gurobi/ ./SCRIPT
```

Here `SCRIPT` is one of the scripts in the folder `neural_network_reachability/examples/tests/` that were created by `make` (e.g., `B_ACC`).

To run all experiments, the script `run_benchmarks.sh` can be used, but the paths to Gurobi need to be adapted. (The last command in the script creates a table, but this is not needed.)



## Plots

For plotting the results after running the experiments, go to the folder `Plots` and run the script `plot.sh`.
