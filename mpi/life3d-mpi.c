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

int lower_bound(int rank, int local_size) {
    return rank*local_size;
}

int upper_bound(int rank, int local_size) {
    return (rank + 1) * local_size;
}

char ***alocate_grid(long long N, int number_of_nodes) {
    int x, y;
    char ***grid;

    int global_N = N * number_of_nodes;
    grid = (char ***) malloc(N * sizeof(char **));
    if (grid == NULL) {
        fprintf(stderr, "Failed to allocate memory for grid\n");
        exit(1);
    }

    #pragma omp parallel for private(x, y) shared(grid)
    for (x = 0; x < N; x++) {
        grid[x] = (char **) malloc(global_N * sizeof(char *));
        if (grid[x] == NULL) {
            fprintf(stderr, "Failed to allocate memory for grid\n");
            exit(1);
        }

        grid[x][0] = (char *) calloc(global_N * global_N, sizeof(char));
        if (grid[x][0] == NULL) {
            fprintf(stderr, "Failed to allocate memory for grid\n");
            exit(1);
        }

        for (y = 1; y < global_N; y++) {
            grid[x][y] = grid[x][0] + y * global_N;
        }
    }

    return grid;
}

char **allocate_frame(int global_N) {
    char **frame = (char **)malloc(global_N * sizeof(char *));
    if (frame == NULL) {
        fprintf(stderr, "Failed to allocate memory for frame\n");
        exit(1);
    }

    char *data = (char *)calloc(global_N * global_N, sizeof(char));
    if (data == NULL) {
        fprintf(stderr, "Failed to allocate memory for data\n");
        exit(1);
    }

    for (int i = 0; i < global_N; i++) {
        frame[i] = data + i * global_N;
    }

    return frame;
}

char ***gen_initial_grid(long long N, float density, int input_seed, int rank, int number_of_nodes, int count_for_first_generation[]) {
    int x, y, z;
    char ***grid;

    int global_N = N * number_of_nodes;
    grid = (char ***) malloc(N * sizeof(char **));
    if (grid == NULL) {
        fprintf(stderr, "Failed to allocate memory for grid\n");
        exit(1);
    }

    #pragma omp parallel for private(x, y, z) shared(grid, density, input_seed)
    for (x = 0; x < N; x++) {
        grid[x] = (char **) malloc(global_N * sizeof(char *));
        if (grid[x] == NULL) {
            fprintf(stderr, "Failed to allocate memory for grid\n");
            exit(1);
        }

        grid[x][0] = (char *) calloc(global_N * global_N, sizeof(char));
        if (grid[x][0] == NULL) {
            fprintf(stderr, "Failed to allocate memory for grid\n");
            exit(1);
        }

        for (y = 1; y < global_N; y++) {
            grid[x][y] = grid[x][0] + y * global_N;
        }
    }

    int species_counts[N_SPECIES + 1] = {0};

    init_r4uni(input_seed);
    // THREAD UNSAFE
    for (x = 0; x < global_N; x++)
        for (y = 0; y < global_N; y++)
            for (z = 0; z < global_N; z++)
                if (r4_uni() < density) {

                    /**
                    convert global X into local X
                    example:
                        - 4 nodes
                        - this process rank = 1
                        - global size = 100
                        we wil only save values from index [25, 49]
                        global X = 36
                        local X = 36 - (25 * 1) = 11
                    */
                    int local_X = x - (N * rank);
                    if (x < upper_bound(rank, N) && x >= lower_bound(rank, N))
                    {
                        grid[local_X][y][z] = (char)(r4_uni() * N_SPECIES) + 1;
                        count_for_first_generation[grid[local_X][y][z]]++;
                    }
                    else
                    {
                        // we call it anyway so that the final seed is the same
                        r4_uni();
                    }
                }
    
    return grid;
}

void free_frame(char **frame, int global_N) {
    if (frame == NULL)
        return;
    free(frame[0]);
    free(frame);
}

void free_grid(char ***grid, int local_N) {
    if (grid == NULL)
        return;

    for (int x = 0; x < local_N; x++) {
        if (grid[x] == NULL)
            continue;
        free(grid[x][0]);
        free(grid[x]);
    }
    free(grid);
}

int count_neighbors(char ***grid, int local_N, int x, int y, int z, char **previous_adjacent_frame, char **next_adjacent_frame, int global_N) {
    int count = 0;
    int dx, dy, dz;

    for (dx = -1; dx <= 1; dx++) {
        for (dy = -1; dy <= 1; dy++) {
            for (dz = -1; dz <= 1; dz++) {
                if (dx == 0 && dy == 0 && dz == 0)
                    continue;
                int nx = (x + dx + local_N) % local_N;
                int ny = (y + dy + global_N) % global_N;
                int nz = (z + dz + global_N) % global_N;
                if (dx == 1 && x == local_N - 1)
                {
                    if (next_adjacent_frame[ny][nz] != 0)
                        count++;
                }
                else if (dx == -1 && x == 0)
                {
                    if (previous_adjacent_frame[ny][nz] != 0)
                        count++;
                }
                else
                {
                    if (grid[nx][ny][nz] != 0)
                    count++;
                }
            }
        }
    }

    return count;
}

void evolve_cell(char ***grid, char ***next_grid, int local_N, int x, int y, int z, int generation, int count_per_generation[], char **previous_adjacent_frame, char **next_adjacent_frame, int global_N) {
    int species = grid[x][y][z];
    int neighbor_count = count_neighbors(grid, local_N, x, y, z, previous_adjacent_frame, next_adjacent_frame, global_N);

    if (species == 0) {  // Empty cell
        if (neighbor_count >= 7 && neighbor_count <= 10) {
            int species_counts[N_SPECIES + 1] = {0};  // Initialize to 0
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dz = -1; dz <= 1; dz++) {
                        if (dx == 0 && dy == 0 && dz == 0)
                            continue;
                        int nx = (x + dx + local_N) % local_N;
                        int ny = (y + dy + global_N) % global_N;
                        int nz = (z + dz + global_N) % global_N;
                        if (dx == 1 && x == local_N - 1)
                        {
                            int neighbor_species = next_adjacent_frame[ny][nz];
                            species_counts[neighbor_species]++;
                        }
                        else if (dx == -1 && x == 0)
                        {
                            int neighbor_species = previous_adjacent_frame[ny][nz];
                            species_counts[neighbor_species]++;
                        }
                        else
                        {
                            int neighbor_species = grid[nx][ny][nz];
                            species_counts[neighbor_species]++;
                        }
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

void send_frames_and_receive(char ***grid, char **previous_adjacent_frame, char **next_adjacent_frame, int local_N, int rank, int number_of_nodes) {
    MPI_Request request;
    int nextRank = (rank + 1) % number_of_nodes;
    int previousRank = (rank - 1 + number_of_nodes) % number_of_nodes; 
    int global_N = local_N * number_of_nodes;
    MPI_Isend(&(grid[local_N - 1][0][0]), global_N * global_N, MPI_CHAR, nextRank, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&(grid[0][0][0]), global_N * global_N, MPI_CHAR, previousRank, 0, MPI_COMM_WORLD, &request);
    MPI_Recv(&(previous_adjacent_frame[0][0]), global_N * global_N, MPI_CHAR, previousRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(next_adjacent_frame[0][0]), global_N * global_N, MPI_CHAR, nextRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
}

void printFrame(char **frame, int global_N) {
    for (int y = 0; y < global_N; y++) {
        for (int z = 0; z < global_N; z++) {
            fprintf(stderr, "%d ", frame[y][z]);
        }
        fprintf(stderr, "\n");
    }
}

void simulation(char ***grid, int local_N, int generations, int rank, int number_of_nodes) {
    char ***next_grid = gen_initial_grid(local_N, 0, seed, rank, number_of_nodes, NULL); // Temporary grid for next generation 
    int global_N = local_N * number_of_nodes;
    char **previous_adjacent_frame;
    char **next_adjacent_frame;

    for (int gen = 1; gen <= generations; gen++) {
        if (number_of_nodes > 1) {
            if (gen == 1) {
                previous_adjacent_frame = allocate_frame(global_N);
                next_adjacent_frame = allocate_frame(global_N);
            }
            send_frames_and_receive(grid, previous_adjacent_frame, next_adjacent_frame, local_N, rank, number_of_nodes);
        } else {
            previous_adjacent_frame = grid[local_N - 1];
            next_adjacent_frame = grid[0];
        }

        // Evolve each cell
        int global_count_per_generation[N_SPECIES + 1] = {0}; // Initialize a global array for reduction
        int root_count_per_generation[N_SPECIES + 1] = {0}; 
        
        #pragma omp parallel shared(grid, next_grid, local_N, gen) 
        {
            int private_count_per_generation[N_SPECIES + 1] = {0}; // Declare private array for each thread

            #pragma omp for collapse(3)
            for (int x = 0; x < local_N; x++) {
                for (int y = 0; y < global_N; y++) {
                    for (int z = 0; z < global_N; z++) {
                        evolve_cell(grid, next_grid, local_N, x, y, z, gen, private_count_per_generation, previous_adjacent_frame, next_adjacent_frame, global_N);
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
            for (int y = 0; y < global_N; y++) {
                for (int z = 0; z < global_N; z++) {
                    grid[x][y][z] = next_grid[x][y][z];
                }
            }
        }

        MPI_Reduce(global_count_per_generation, root_count_per_generation, N_SPECIES + 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            // Update max_species_count_list and max_generation using global_count_per_generation
            #pragma omp parallel for
            for (int i = 1; i <= N_SPECIES; i++) {
                #pragma omp critical
                {
                    if (root_count_per_generation[i] > max_species_count_list[i]) {
                        max_species_count_list[i] = root_count_per_generation[i];
                        max_generation[i] = gen;
                    }
                }
            }
        }

    }

    if (number_of_nodes > 1)
    {
        free_frame(previous_adjacent_frame, global_N);
        free_frame(next_adjacent_frame, global_N);
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
    int count_for_first_generation[N_SPECIES + 1] = {0}; 
    int root_count_for_first_generation[N_SPECIES + 1] = {0};
    grid = gen_initial_grid(local_N, density, seed, rank, size, count_for_first_generation); // Generate initial grid for each process

    MPI_Reduce(count_for_first_generation, root_count_for_first_generation, N_SPECIES + 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        // Update max_species_count_list and max_generation using global_count_per_generation
        #pragma omp parallel for
        for (int i = 1; i <= N_SPECIES; i++) {
            #pragma omp critical
            {
                if (root_count_for_first_generation[i] > max_species_count_list[i]) {
                    max_species_count_list[i] = root_count_for_first_generation[i];
                    max_generation[i] = 0;
                }
            }
        }
    }
    
    exec_time = -MPI_Wtime();

    // Simulate game of life
    simulation(grid, local_N, generations, rank, size);

    exec_time += MPI_Wtime();

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
