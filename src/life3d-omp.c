#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N_SPECIES 9

unsigned int seed;

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

        for (y = 1; y < N; y++)
            grid[x][y] = grid[x][0] + y * N;
    }

    init_r4uni(input_seed);
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            for (z = 0; z < N; z++)
                if (r4_uni() < density)
                    grid[x][y][z] = (int)(r4_uni() * N_SPECIES) + 1;

    return grid;
}

void free_grid(char ***grid, long long N) {
    if (grid == NULL)
        return;

    for (int x = 0; x < N; x++) {
        if (grid[x] == NULL)
            continue;
        free(grid[x][0]);
        free(grid[x]);
    }
    free(grid);
}

int count_neighbors(char ***grid, long long N, int x, int y, int z, int species) {
    int count = 0;
    int dx, dy, dz;

    for (dx = -1; dx <= 1; dx++) {
        for (dy = -1; dy <= 1; dy++) {
            for (dz = -1; dz <= 1; dz++) {
                if (dx == 0 && dy == 0 && dz == 0)
                    continue;
                int nx = (x + dx + N) % N;
                int ny = (y + dy + N) % N;
                int nz = (z + dz + N) % N;
                if (grid[nx][ny][nz] == species)
                    count++;
            }
        }
    }

    return count;
}

void evolve_cell(char ***grid, char ***next_grid, long long N, int x, int y, int z) {
    int species = grid[x][y][z];
    int neighbor_count = count_neighbors(grid, N, x, y, z, species);
    
    if (species == 0) {  // Empty cell
        if (neighbor_count >= 7 && neighbor_count <= 10) {
            // Reproduction
            int species_counts[N_SPECIES + 1] = {0};  // Initialize to 0
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dz = -1; dz <= 1; dz++) {
                        if (dx == 0 && dy == 0 && dz == 0)
                            continue;
                        int nx = (x + dx + N) % N;
                        int ny = (y + dy + N) % N;
                        int nz = (z + dz + N) % N;
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
            // Underpopulation or Overcrowding
            next_grid[x][y][z] = 0;  // Die
        } else {
            // Survive
            next_grid[x][y][z] = species;
        }
    }
}

void simulation(char ***grid, long long N, int generations) {
    char ***next_grid = gen_initial_grid(N, 0.0, 0);  // Temporary grid for next generation
    
    for (int gen = 0; gen < generations; gen++) {
        // Evolve each cell
        for (int x = 0; x < N; x++) {
            for (int y = 0; y < N; y++) {
                for (int z = 0; z < N; z++) {
                    evolve_cell(grid, next_grid, N, x, y, z);
                }
            }
        }
        // Copy next_grid back to grid for next generation
        for (int x = 0; x < N; x++) {
            for (int y = 0; y < N; y++) {
                for (int z = 0; z < N; z++) {
                    grid[x][y][z] = next_grid[x][y][z];
                }
            }
        }
    }
    free_grid(next_grid, N);
}

void print_result(char ***grid, long long N) {
    // Count population for each species over all generations
    int max_species_count[N_SPECIES + 1] = {0};  // Initialize to 0
    int max_generation[N_SPECIES + 1] = {0};     // Initialize to 0
    
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                int species = grid[x][y][z];
                if (species != 0) {  // Ignore empty cells
                    max_species_count[species]++;
                    if (max_generation[species] == 0) {
                        max_generation[species] = 1;  // Initialize to first generation
                    }
                }
            }
        }
    }
    
    // Print population statistics for each species
    for (int i = 1; i <= N_SPECIES; i++) {
        printf("%d %d %d\n", i, max_species_count[i], max_generation[i]);
    }
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <generations> <cube_size> <density> <seed>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    int generations = atoi(argv[1]);
    long long N = atoll(argv[2]); // Use long long for cube size
    float density = atof(argv[3]);
    int seed = atoi(argv[4]);

    // Generate initial grid
    char ***grid = gen_initial_grid(N, density, seed);

    // Simulate game of life
    simulation(grid, N, generations);

    // Output results
    print_result(grid, N);

    // Free memory
    free_grid(grid, N);

    return 0;
}
