#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ALPHA 0.01
#define DX 0.02
#define DY 0.02
#define DT 0.0005
#define T 1500
#define BMP_HEADER_SIZE 54

void write_bmp_header(FILE *file, int width, int height) {
    unsigned char header[BMP_HEADER_SIZE] = {0};
    int file_size = BMP_HEADER_SIZE + 3 * width * height;
    header[0] = 'B'; header[1] = 'M';
    header[2] = file_size & 0xFF; header[3] = (file_size >> 8) & 0xFF;
    header[4] = (file_size >> 16) & 0xFF; header[5] = (file_size >> 24) & 0xFF;
    header[10] = BMP_HEADER_SIZE;
    header[14] = 40;
    header[18] = width & 0xFF; header[19] = (width >> 8) & 0xFF;
    header[20] = (width >> 16) & 0xFF; header[21] = (width >> 24) & 0xFF;
    header[22] = height & 0xFF; header[23] = (height >> 8) & 0xFF;
    header[24] = (height >> 16) & 0xFF; header[25] = (height >> 24) & 0xFF;
    header[26] = 1; header[28] = 24;
    fwrite(header, 1, BMP_HEADER_SIZE, file);
}

void get_color(double value, unsigned char *r, unsigned char *g, unsigned char *b) {
    if (value >= 500.0) {*r = 255; *g = 0; *b = 0;}
    else if (value >= 100.0) {*r = 255; *g = 128; *b = 0;}
    else if (value >= 50.0) {*r = 171; *g = 71; *b = 188;}
    else if (value >= 25) {*r = 255; *g = 255; *b = 0;}
    else if (value >= 1) {*r = 0; *g = 0; *b = 255;}
    else if (value >= 0.1) {*r = 5; *g = 248; *b = 252;}
    else {*r = 255; *g = 255; *b = 255;}
}

void write_grid(FILE *file, double *grid, int nx, int ny) {
    for (int i = nx - 1; i >= 0; i--) {
        for (int j = 0; j < ny; j++) {
            unsigned char r, g, b;
            get_color(grid[i * ny + j], &r, &g, &b);
            fwrite(&b, 1, 1, file);
            fwrite(&g, 1, 1, file);
            fwrite(&r, 1, 1, file);
        }
        for (int pad = 0; pad < (4 - (ny * 3) % 4) % 4; pad++) {
            fputc(0, file);
        }
    }
}

void initialize_grid(double *grid, int nx, int ny) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            if (i == j || i == (ny - 1 - j))
                grid[i * ny + j] = T;
            else
                grid[i * ny + j] = 0.0;
        }
    }
}

void solve_heat(double *grid, double *new_grid, int local_nx, int ny, double r, int rank, int size) {
    MPI_Status status;
    // Exchange boundaries
    if (rank > 0) {
        MPI_Sendrecv(&grid[1 * ny], ny, MPI_DOUBLE, rank - 1, 0,
                     &grid[0 * ny], ny, MPI_DOUBLE, rank - 1, 0,
                     MPI_COMM_WORLD, &status);
    }
    if (rank < size - 1) {
        MPI_Sendrecv(&grid[local_nx * ny], ny, MPI_DOUBLE, rank + 1, 0,
                     &grid[(local_nx + 1) * ny], ny, MPI_DOUBLE, rank + 1, 0,
                     MPI_COMM_WORLD, &status);
    }

    // Compute
    #pragma omp parallel for collapse(2)
    for (int i = 1; i <= local_nx; i++) {
        for (int j = 1; j < ny - 1; j++) {
            new_grid[i * ny + j] = grid[i * ny + j]
                + r * (grid[(i + 1) * ny + j] + grid[(i - 1) * ny + j] - 2 * grid[i * ny + j])
                + r * (grid[i * ny + j + 1] + grid[i * ny + j - 1] - 2 * grid[i * ny + j]);
        }
    }

    // Dirichlet boundaries
    if (rank == 0) {
        #pragma omp parallel for
        for (int j = 0; j < ny; j++) new_grid[1 * ny + j] = 0.0;
    }
    if (rank == size - 1) {
        #pragma omp parallel for
        for (int j = 0; j < ny; j++) new_grid[local_nx * ny + j] = 0.0;
    }
    #pragma omp parallel for
    for (int i = 1; i <= local_nx; i++) {
        new_grid[i * ny + 0] = 0.0;
        new_grid[i * ny + (ny - 1)] = 0.0;
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == 0)
            fprintf(stderr, "Usage: %s <grid_size> <steps> <num_threads>\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int nx = atoi(argv[1]);
    int ny = nx;
    int steps = atoi(argv[2]);
    int num_threads = atoi(argv[3]);
    omp_set_num_threads(num_threads);

    double r = ALPHA * DT / (DX * DY);
    if (r > 0.25 && rank == 0)
        printf("WARNING: Scheme may be unstable (r = %.3f > 0.25)\n", r);

    // Compute local size and displacement
    int base = nx / size, rem = nx % size;
    int local_nx = base + (rank < rem ? 1 : 0);
    int start_row = rank * base + (rank < rem ? rank : rem);

    double *grid = calloc((local_nx + 2) * ny, sizeof(double));
    double *new_grid = calloc((local_nx + 2) * ny, sizeof(double));

    double *global_grid = NULL;
    int *counts = NULL, *displs = NULL;

    if (rank == 0) {
        global_grid = malloc(nx * ny * sizeof(double));
        initialize_grid(global_grid, nx, ny);

        counts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        for (int i = 0, offset = 0; i < size; i++) {
            int rows = base + (i < rem ? 1 : 0);
            counts[i] = rows * ny;
            displs[i] = offset;
            offset += counts[i];
        }
    }

    MPI_Scatterv(global_grid, counts, displs, MPI_DOUBLE,
                 &grid[1 * ny], local_nx * ny, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    double t0 = MPI_Wtime();
    for (int s = 0; s < steps; s++) {
        solve_heat(grid, new_grid, local_nx, ny, r, rank, size);
        double *tmp = grid; grid = new_grid; new_grid = tmp;
    }
    double t1 = MPI_Wtime();

    if (rank == 0) printf("Hybrid Execution Time (MPI+OpenMP): %.6f seconds\n", t1 - t0);

    MPI_Gatherv(&grid[1 * ny], local_nx * ny, MPI_DOUBLE,
                global_grid, counts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        char filename[256];
        snprintf(filename, sizeof(filename), "output_hybrid_%d_%d_mpi%d_omp%d.bmp", nx, steps, size, num_threads);
        FILE *f = fopen(filename, "wb");
        if (!f) {
            perror("Erreur d'ouverture du fichier");
        } else {
            write_bmp_header(f, ny, nx);
            write_grid(f, global_grid, nx, ny);
            fclose(f);
            printf("Wrote %s\n", filename);
        }
        free(global_grid);
        free(counts);
        free(displs);
    }

    free(grid);
    free(new_grid);
    MPI_Finalize();
    return 0;
}
