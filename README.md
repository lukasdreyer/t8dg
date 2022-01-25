### Introduction

t8dg (spoken as "tet-de-ge") is a C/C++ software to use a discontinuous galerkin method on adaptive meshes.
t8dg uses t8code as a mesh managing library. 
It is licensed under the GNU General Public License 3.0 or later. Copyright (c) 2020 the developers.


<table>
    <tr>
        <td><img src="https://github.com/lukasdreyer/t8dg/blob/main/doc/pictures/adapt_2D.png?raw=true" height="300" /></td> 
        <td><img src="https://github.com/lukasdreyer/t8dg/blob/main/doc/pictures/cylinder_cut83.png?raw=true" height="300" /></td>
    </tr>
</table>

### Setup

We provide a short guide to install t8dg. 

For a more detailed description, please see the [Installation guide](https://github.com/lukasdreyer/t8dg/blob/main/INSTALL).

#### Requirements

- [libsc](https://github.com/cburstedde/libsc) (Included in t8dg's git repository)
- [p4est](https://github.com/cburstedde/p4est) (Included in t8dg's git repository)
- [t8code](https://github.com/holke/t8code) (Included in t8dg's git repository)
- automake
- libtool
- make

#### Steps
To setup the project perform the following steps
  
    1.) If you cloned from github, initialize and download the git submodules
       t8code, p4est and sc.
      - git submodule init
      - git submodule update      
    2.) Call the bootstrap script in the source directory
      - ./bootstrap        
    3.) Goto your installation folder and call configure and make
      - cd /path/to/install
      - /path/to/source/configure [OPTIONS]
      - make 
      - make check
      - make install

To see a list of possible configure options, call
 
 ./configure -h

Most commonly used for t8dg are

  --enable-mpi    (enables MPI parallelization)
  
  --enable-debug  (enables debugging mode - massively reduces performance)
  
#### Example configurations

For a parallel release mode with local installation path `$HOME/t8dg_install`:

`configure --enable-mpi CFLAGS=-O3 CXXFLAGS=-O3 --prefix=$HOME/t8dg_install`

For a debugging mode with static linkage (makes using gdb and valgrind more comfortable):

`configure --enable-mpi --enable-debug --enable-static --disable-shared CFLAGS="-Wall -O0 -g" CXXFLAGS="-Wall -O0 -g"`
  
### Getting started
      1. configure and make as above
      2. Go to exec folder
        - cd /path/to/install
        - cd exec
      3. Call t8dg with -h flag to show possible options 
        - ./t8dg -h
      4. Call t8dg with appropriate options, for example with 
        - minimum refinement level 3, maximum refinement level 5 (-l3 -r5)
        - spatial and time integration order 4 (-L4 -o4)
        - CFL number = 0.05 (-C0.05) (can be chosen higher for smaller integration order)
        - square mesh (-m3)
        - and circular step function as initial condition (-i6):
        - ./t8dg -l3 -r5 -L4 -o4 -C0.05 -m3 -i6
      5. If t8dg is configured with MPI, call it with mpirun to execute in parallel:
        - mpirun -n 4 ./t8dg -l3 -r5 -L4 -o4 -C0.05 -m3 -i6
      6. t8dg produces the vtu files 't8dg_advection_*.{p}vtu'.
         To investigate, you can open the pvtu files in paraview.
         We suggest coloring for 'Num. Solution' and enabling 'Surface with edges'.

  ### Publications
  
  An (incomplete) list of publications related to t8dg:
    
  [1] Lukas Dreyer, The local discontinuous galerkin method for the advection-diffusion equation on adaptive meshes, Master thesis at University of Bonn, 2021,
      [Full text available](https://elib.dlr.de/143969/)
      
  
  ### Citing t8dg
  
  If you use t8dg in any of your publications, please cite the [github repository](https://github.com/lukasdreyer/t8dg) and [1].
