#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define BMP_HEADER_SIZE 54
#define ALPHA 0.01 // Thermal diffusivity
#define L 0.2      // Length (m) of the square domain
#define DX 0.02    // grid spacing in x-direction
#define DY 0.02    // grid spacing in y-direction
#define DT 0.0005  // Time step
#define T 1500     // Temperature on ºk of the heat source

// Function to print the grid (optional, for debugging or visualization)
void print_grid(double *grid, int nx, int ny) {
  int i, j;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      printf("%.2f ", grid[i * ny + j]);
    }
    printf("\n");
  }
  printf("\n");
}
// Function to initialize the grid
void initialize_grid(double *grid, int nx, int ny, int temp_source) {
  int i, j;
#pragma omp parallel for private(i, j) collapse(2)
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      if (i == j)
        grid[i * ny + j] = 1500.0;
      else if (i == nx - 1 - j)
        grid[i * ny + j] = 1500.0;
      else
        grid[i * ny + j] = 0.0;
    }
  }
}

void solve_heat_equation_hybrid(double *grid, double *new_grid, int steps,
                                double r, int nx, int ny, int local_nx,
                                int start_i, int rank, int size) {
  int step, i, j;
  double *temp;
  MPI_Status status;

  for (step = 0; step < steps; step++) {
    // 1) Exchange border rows: send first real row to rank-1, recv into row 0;
    //    send last real row to rank+1, recv into row local_nx+1.
    if (rank > 0) {
      MPI_Sendrecv(&grid[1 * ny], ny, MPI_DOUBLE, rank - 1, 0, &grid[0 * ny],
                   ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
    }
    if (rank < size - 1) {
      MPI_Sendrecv(&grid[local_nx * ny], ny, MPI_DOUBLE, rank + 1, 0,
                   &grid[(local_nx + 1) * ny], ny, MPI_DOUBLE, rank + 1, 0,
                   MPI_COMM_WORLD, &status);
    }

// Compute the new grid
#pragma omp parallel for collapse(2) private(i, j)
    for (i = 1; i <= local_nx; i++) {
      for (j = 1; j < ny - 1; j++) {
        new_grid[i * ny + j] =
            grid[i * ny + j] +
            r * (grid[(i + 1) * ny + j] + grid[(i - 1) * ny + j] -
                 2 * grid[i * ny + j]) +
            r * (grid[i * ny + j + 1] + grid[i * ny + j - 1] -
                 2 * grid[i * ny + j]);
      }
    }

    // Apply boundary conditions (Dirichlet: u=0 on boundaries)
    //    Top boundary only on rank 0’s first real row (i==1 corresponds to
    //    global i==0)
    if (rank == 0) {
#pragma omp parallel for
      for (j = 0; j < ny; j++) {
        new_grid[1 * ny + j] = 0.0;
      }
    }
    //    - Bottom boundary only on rank size-1’s last real row
    if (rank == size - 1) {
#pragma omp parallel for
      for (j = 0; j < ny; j++) {
        new_grid[local_nx * ny + j] = 0.0;
      }
    }
//    - Left and right boundaries on every rank’s subdomain
#pragma omp parallel for private(i)
    for (i = 1; i <= local_nx; i++) {
      new_grid[i * ny + 0] = 0.0;
      new_grid[i * ny + (ny - 1)] = 0.0;
    }

    // Swap the grids
    temp = grid;
    grid = new_grid;
    new_grid = temp;
  }
}
// Function to write BMP file header
void write_bmp_header(FILE *file, int width, int height) {
  unsigned char header[BMP_HEADER_SIZE] = {0};

  int file_size = BMP_HEADER_SIZE + 3 * width * height;
  header[0] = 'B';
  header[1] = 'M';
  header[2] = file_size & 0xFF;
  header[3] = (file_size >> 8) & 0xFF;
  header[4] = (file_size >> 16) & 0xFF;
  header[5] = (file_size >> 24) & 0xFF;
  header[10] = BMP_HEADER_SIZE;

  header[14] = 40; // Info header size
  header[18] = width & 0xFF;
  header[19] = (width >> 8) & 0xFF;
  header[20] = (width >> 16) & 0xFF;
  header[21] = (width >> 24) & 0xFF;
  header[22] = height & 0xFF;
  header[23] = (height >> 8) & 0xFF;
  header[24] = (height >> 16) & 0xFF;
  header[25] = (height >> 24) & 0xFF;
  header[26] = 1;  // Planes
  header[28] = 24; // Bits per pixel

  fwrite(header, 1, BMP_HEADER_SIZE, file);
}

void get_color(double value, unsigned char *r, unsigned char *g,
               unsigned char *b) {

  if (value >= 500.0) {
    *r = 255;
    *g = 0;
    *b = 0; // Red
  } else if (value >= 100.0) {
    *r = 255;
    *g = 128;
    *b = 0; // Orange
  } else if (value >= 50.0) {
    *r = 171;
    *g = 71;
    *b = 188; // Lilac
  } else if (value >= 25) {
    *r = 255;
    *g = 255;
    *b = 0; // Yellow
  } else if (value >= 1) {
    *r = 0;
    *g = 0;
    *b = 255; // Blue
  } else if (value >= 0.1) {
    *r = 5;
    *g = 248;
    *b = 252; // Cyan
  } else {
    *r = 255;
    *g = 255;
    *b = 255; // white
  }
}
// Function to write the grid matrix into the file
void write_grid(FILE *file, double *grid, int nx, int ny) {
  int i, j, padding;
  // Write pixel data to BMP file
  for (i = nx - 1; i >= 0; i--) { // BMP format stores pixels bottom-to-top
    for (j = 0; j < ny; j++) {
      unsigned char r, g, b;
      get_color(grid[i * ny + j], &r, &g, &b);
      fwrite(&b, 1, 1, file); // Write blue channel
      fwrite(&g, 1, 1, file); // Write green channel
      fwrite(&r, 1, 1, file); // Write red channel
    }
    // Row padding for 4-byte alignment (if necessary)
    for (padding = 0; padding < (4 - (nx * 3) % 4) % 4; padding++) {
      fputc(0, file);
    }
  }
}

int main(int argc, char *argv[]) {
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int nx, ny, steps;
  double r;
  char filename[256];
  if (rank == 0) {
    if (argc != 4) {
      fprintf(stderr, "Usage: %s size steps name_output_file.bmp\n", argv[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    nx = atoi(argv[1]);
    ny = nx;
    steps = atoi(argv[2]);
    snprintf(filename, sizeof(filename), "%s", argv[3]);
    r = ALPHA * DT / (DX * DY);
  }

  // send the parameters to all ranks
  MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&r, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // compute each rank grid
  int base = nx / size;
  int rem = nx % size;
  int local_nx = base + (rank < rem ? 1 : 0);
  int start_i = rank * base + (rank < rem ? rank : rem);

  // allocate local arrays (with 2 ghost‐rows)
  double *local_grid = calloc((local_nx + 2) * ny, sizeof(double));
  double *local_new_grid = calloc((local_nx + 2) * ny, sizeof(double));
  if (!local_grid || !local_new_grid) {
    fprintf(stderr, "Rank %d: allocation failed\n", rank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // prepare scatter/gather counts & displacements
  int *counts = malloc(size * sizeof(int));
  int *displs = malloc(size * sizeof(int));
  int rnk;
  for (rnk = 0; rnk < size; rnk++) {
    int ln = base + (rnk < rem ? 1 : 0);
    counts[rnk] = ln * ny;
    displs[rnk] = ((rnk * base) + (rnk < rem ? rnk : rem)) * ny;
  }

  // rank 0 initializes the full grid
  double *full_grid = NULL;
  if (rank == 0) {
    full_grid = malloc(nx * ny * sizeof(double));
    initialize_grid(full_grid, nx, ny, T);
  }
  // scatter into grid[1..local_nx] (first row & local_nx+1 are ghosts)
  MPI_Scatterv(full_grid, counts, displs, MPI_DOUBLE, &local_grid[1 * ny],
               local_nx * ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    free(full_grid);
    full_grid = NULL;
  }

  // sync and start timer
  MPI_Barrier(MPI_COMM_WORLD);
  double start_time = MPI_Wtime();

  solve_heat_equation_hybrid(local_grid, local_new_grid, steps, r, nx, ny,
                             local_nx, start_i, rank, size);

  // wait for all ranks to finish
  MPI_Barrier(MPI_COMM_WORLD);
  double end_time = MPI_Wtime();

  // gather the interior rows back to full_grid on rank 0
  if (rank == 0) {
    full_grid = malloc(nx * ny * sizeof(double));
  }
  MPI_Gatherv(&local_grid[1 * ny], local_nx * ny, MPI_DOUBLE, full_grid, counts,
              displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // rank 0 writes out .bmp and prints timing
  if (rank == 0) {
    FILE *f = fopen(filename, "wb");
    if (!f) {
      fprintf(stderr, "Error opening output file \"%s\"\n", filename);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    write_bmp_header(f, nx, ny);
    write_grid(f, full_grid, nx, ny);
    fclose(f);
    printf("Execution Time = %f seconds\n", end_time - start_time);
    free(full_grid);
  }

  free(local_grid);
  free(local_new_grid);
  free(counts);
  free(displs);

  MPI_Finalize();
  return 0;
}
