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

To run the program, use the following command:

```bash
./heat_serial <matrix_size> <steps> <output_file.bmp>
```

Replace `<matrix_size>` with the size of the square matrix (e.g., 100, 1000), `<steps>` with the number of time steps, and `<output_file.bmp>` with the desired output BMP file name.

# TODO

- [ ] Report
- [ ] Analyze the serial code and detect the code susceptible to be parallelized.
- [ ] Justify the strategy used to parallelize the serial code. You should point
      out how has solved the data dependencies.
- [ ] Analyse if it is possible to parallelize the output and justify the answer.
- [ ] If you have tested different strategies, you can describe briefly each of
      them and justify the reason of your choice. It is not needed to test for all
      the test-bed of data.
- [ ] Analyze scalability and speedup in relation to the number of threads and
      size problem. Consider on your discussion the effects of the dependences
      between data, in case it was needed, in your solution.

## Serial results

| Matrix Size \ Steps | 100 | 1000 | 10000 | 100000 |
| ------------------- | --- | ---- | ----- | ------ |
| 100x100             |     |      |       |        |
| 1000x1000           |     |      |       |        |
| 2000x2000           |     |      |       |        |

## OpenMP results

| Matrix Size \ Steps | 100 | 1000 | 10000 | 100000 |
| ------------------- | --- | ---- | ----- | ------ |
| 100x100             |     |      |       |        |
| 1000x1000           |     |      |       |        |
| 2000x2000           |     |      |       |        |
