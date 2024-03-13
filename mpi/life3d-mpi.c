#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#define N_SPECIES 9

unsigned int seed;

int max_species_count_list[N_SPECIES + 1] = {0};  // Initialize to 0
int max_generation[N_SPECIES + 1] = {0};     // Initialize to 0

void init_r4uni(int input_seed) {
    seed = input_seed + 987654321;
}

float r4_uni() {
    int seed_in = seed;

    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);

    return 0.5 + 0.2328306e-09 * (seed_in + (int) seed);
}

char ***gen_initial_grid(long long N, float density, int input_seed) {
    int x, y, z;
    char ***grid;

    grid = (char ***) malloc(N * sizeof(char **));
    if (grid == NULL) {
        fprintf(stderr, "Failed to allocate memory for grid\n");
        exit(1);
    }

    #pragma omp parallel for private(x, y, z) shared(grid, density, input_seed)
    for (x = 0; x < N; x++) {
        grid[x] = (char **) malloc(N * sizeof(char *));
        if (grid[x] == NULL) {
            fprintf(stderr, "Failed to allocate memory for grid\n");
            exit(1);
        }

        grid[x][0] = (char *) calloc(N * N, sizeof(char));
        if (grid[x][0] == NULL) {
            fprintf(stderr, "Failed to allocate memory for grid\n");
            exit(1);
        }

        for (y = 1; y < N; y++) {
            grid[x][y] = grid[x][0] + y * N;
        }
    }

    int species_counts[N_SPECIES + 1] = {0};

    init_r4uni(input_seed);
    // THREAD UNSAFE
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            for (z = 0; z < N; z++)
                if (r4_uni() < density) {
                    grid[x][y][z] = (char)(r4_uni() * N_SPECIES) + 1;
                    #pragma omp atomic
                    species_counts[grid[x][y][z]]++;
                    if (species_counts[grid[x][y][z]] > max_species_count_list[grid[x][y][z]]) 
                    {   
                        #pragma omp critical
                        {
                            max_species_count_list[grid[x][y][z]] = species_counts[grid[x][y][z]];
                            max_generation[grid[x][y][z]] = 0;
                        }
                    }
                }
    
    return grid;
}

void free_grid(char ***grid, int local_N) {
    if (grid == NULL)
        return;

    #pragma omp parallel for
    for (int x = 0; x < local_N; x++) {
        if (grid[x] == NULL)
            continue;
        free(grid[x][0]);
        free(grid[x]);
    }

    free(grid);
}

int count_neighbors(char ***grid, int local_N, int x, int y, int z) {
    int count = 0;
    int dx, dy, dz;

    for (dx = -1; dx <= 1; dx++) {
        for (dy = -1; dy <= 1; dy++) {
            for (dz = -1; dz <= 1; dz++) {
                if (dx == 0 && dy == 0 && dz == 0)
                    continue;
                int nx = (x + dx + local_N) % local_N;
                int ny = (y + dy + local_N) % local_N;
                int nz = (z + dz + local_N) % local_N;
                if (grid[nx][ny][nz] != 0)
                    count++;
            }
        }
    }

    return count;
}

void evolve_cell(char ***grid, char ***next_grid, int local_N, int x, int y, int z, int generation, int count_per_generation[] ) {
    int species = grid[x][y][z];
    int neighbor_count = count_neighbors(grid, local_N, x, y, z);

    if (species == 0) {  // Empty cell
        if (neighbor_count >= 7 && neighbor_count <= 10) {
            int species_counts[N_SPECIES + 1] = {0};  // Initialize to 0
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dz = -1; dz <= 1; dz++) {
                        if (dx == 0 && dy == 0 && dz == 0)
                            continue;
                        int nx = (x + dx + local_N) % local_N;
                        int ny = (y + dy + local_N) % local_N;
                        int nz = (z + dz + local_N) % local_N;
                        int neighbor_species = grid[nx][ny][nz];
                        species_counts[neighbor_species]++;
                    }
                }
            }
            int max_species_count = 0;
            int dominant_species = 0;
            for (int i = 1; i <= N_SPECIES; i++) {
                if (species_counts[i] > max_species_count) {
                    max_species_count = species_counts[i];
                    dominant_species = i;
                }
            }
            next_grid[x][y][z] = dominant_species;
        } else {
            next_grid[x][y][z] = 0;  // Remain empty
        }
    } else {  // Occupied cell
        if (neighbor_count <= 4 || neighbor_count > 13) {
            next_grid[x][y][z] = 0;  // Die
        } else {
            next_grid[x][y][z] = species;  // Survive
        }
    }

    count_per_generation[next_grid[x][y][z]]++;
}

void simulation(char ***grid, int local_N, int generations) {
    char ***next_grid = gen_initial_grid(local_N, 0.0, 0);  // Temporary grid for next generation

    for (int gen = 1; gen <= generations; gen++) {
        // Evolve each cell
        int global_count_per_generation[N_SPECIES + 1] = {0}; // Initialize a global array for reduction
        
        #pragma omp parallel shared(grid, next_grid, local_N, gen) 
        {
            int private_count_per_generation[N_SPECIES + 1] = {0}; // Declare private array for each thread

            #pragma omp for collapse(3)
            for (int x = 0; x < local_N; x++) {
                for (int y = 0; y < local_N; y++) {
                    for (int z = 0; z < local_N; z++) {
                        evolve_cell(grid, next_grid, local_N, x, y, z, gen, private_count_per_generation);
                    }
                }
            }

            // Accumulate private counts into global array using a critical section
            #pragma omp critical 
            {
                for (int i = 1; i <= N_SPECIES; i++) {
                    global_count_per_generation[i] += private_count_per_generation[i];
                }
            }
        }
        
        // Copy next_grid back to grid for next generation
        #pragma omp parallel for
        for (int x = 0; x < local_N; x++) {
            for (int y = 0; y < local_N; y++) {
                for (int z = 0; z < local_N; z++) {
                    grid[x][y][z] = next_grid[x][y][z];
                }
            }
        }

        // Update max_species_count_list and max_generation using global_count_per_generation
        #pragma omp parallel for
        for (int i = 1; i <= N_SPECIES; i++) {
            #pragma omp critical
            {
                if (global_count_per_generation[i] > max_species_count_list[i]) {
                    max_species_count_list[i] = global_count_per_generation[i];
                    max_generation[i] = gen;
                }
            }
        }
    }

    free_grid(next_grid, local_N);
}

void print_result(int rank, int size) {
    // Print population statistics for each species
    for (int i = 1; i <= N_SPECIES; i++) {
        printf("%d %d %d\n", i, max_species_count_list[i], max_generation[i]);
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    double exec_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 5) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <generations> <cube_size> <density> <seed>\n", argv[0]);
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    int generations = atoi(argv[1]);
    int global_N = atoi(argv[2]); // Total size of the cube
    float density = atof(argv[3]);
    int seed = atoi(argv[4]);

    if (generations <= 0 || global_N <= 0 || density < 0 || density > 1) {
        if (rank == 0) {
            fprintf(stderr, "Invalid input\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    int local_N = global_N / size; // Size of the cube for each process
    char ***grid;
    
    if (rank == 0) {
        grid = gen_initial_grid(local_N, density, seed + rank); // Generate initial grid for each process
    }

    // Scatter initial grid from root to all processes
    // MPI_Scatter(...);

    //MPI_Scatter(grid[0][0], local_N * local_N * local_N, MPI_CHAR, grid[0][0], local_N * local_N * local_N, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Bcast(grid[0][0], local_N * local_N * local_N, MPI_CHAR, 0, MPI_COMM_WORLD);

    exec_time = -MPI_Wtime();

    // Simulate game of life
    simulation(grid, local_N, generations);

    exec_time += MPI_Wtime();

    MPI_Reduce(MPI_IN_PLACE, max_species_count_list, N_SPECIES + 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, max_generation, N_SPECIES + 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    // Gather max_species_count_list and max_generation from all processes to root
    // MPI_Gather(...);
    //MPI_Gather(max_species_count_list, N_SPECIES + 1, MPI_INT, max_species_count_list, N_SPECIES + 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        fprintf(stderr, "%.1fs\n", exec_time);
        print_result(rank, size);
    }

    if (rank == 0) {
        free_grid(grid, local_N);
    }

    MPI_Finalize();

    return 0;
}
