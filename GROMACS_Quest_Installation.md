Using the latest gromacs version will always be the fastest but custom optmized install is a little tricky. 
Here is an example installation procedure to help.

## Installation directions of GROMACS 2020.3 on Quest with MPI and CUDA enabled
We will be installing GROMACS in the local user directory `/home/<user>/.local`. 
`mkdir /home/<user>/.local` directory if you don't have already have a local directory to install custom software. If you already have such a local directory, you can use that.
`<user>` is your username.

Download and unzip GROMACS 2020.3 from https://manual.gromacs.org/documentation/

Load essential modules.
```bash
module purge all
module load cmake/3.12.0
module load cuda/cuda-9.2
```
You are also put the above `cmake` and `cuda` load commands in your `~/.bash_profile` so that these module load automatically when you ssh to the quest server.

```
cd <path>/<to>/gromacs-2020.3
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$HOME/.local/gromacs -DGMX_FFT_LIBRARY=fftw3 -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which c++) -DGMX_GPU=on -DGMX_MPI=on -DPYTHON_EXECUTABLE=/usr/bin/python3
make -j8
make install
```

Finally add the following line in your `~/.bash_profile`
`PATH="$HOME/.local/gromacs/bin:$PATH"`


### A quick benchmark using this GROMACS version
An atomistic simulation consisting 200 PAs in a 12x12x12 nm box for 100 ns completes in 
- 7 days using 1 GPU  | command: gmx mdrun -deffnm md -v -cpt 5 -nb gpu -pme gpu -update gpu
- 4 days using 4 GPUs | command: mpirun -np 4 gmx_mpi mdrun -deffnm md -v -cpt 5 -nb gpu -pme gpu -npme 1
