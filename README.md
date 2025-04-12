# Heat Diffusion Equation

This project aims to solve the Heat Diffusion Equation problem by means of OpenMP.
Specifically, the main aim is to parallelize the given problem and analyse the
scalability of the proposed solutions and the effects of design decisions according
to the characteristics of the problem.

## Build and Run Instructions

To build the project, use the following command:

```bash
make
```
It will compile both parallel and serial implementations that can be also compiled using `make parallel` or `make serial`. 

To run the program, use the following command:

```bash
./heat_serial <matrix_size> <steps> <output_file.bmp>
```

Replace `<matrix_size>` with the size of the square matrix (e.g., 100, 1000), `<steps>` with the number of time steps, and `<output_file.bmp>` with the desired output BMP file name.
