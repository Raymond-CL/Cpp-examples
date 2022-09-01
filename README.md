## C++ Examples

some C-plus-plus examples

---

### Conway's Game of Life

Cellular Automata Simuation
In memory of John Horton Conway's passing on 11th April, 2020

Rules:
- live cell with 2 or 3 live neighbours lives on, else dies.
- dead cell with 3 live neighbours lives.

Comments:
- run on both unix or windows OS
- randomized initialization
- flashing the display screen is kinda slow for different OS

### gsl

an example to test the GNU Scientific Library (GSL) 
- to install GSL on Ubuntu WSL, follow 
```
sudo apt install libgsl-dev
```
- to link the GSL library, use the `-lgsl -lm` compiler option
- the example calculates the volume of a hyper-sphere with given radius and dimension
- the example utilizes the `plain`, `miser` and `vegas` Monte-Carlo integration
- the example also utilizes some GSL math functions, such as the Gamma function

