#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char **argv) {
    int i, j, k, l, rc, rank, p, choice, N, n, *temp;
    int *A, *B, **C, *elements_C, **D, *elements_D, **CplusD, *elements_CplusD;
    int *CprodB, local_AprodB, AprodB, source_p, dest_p, localD_index, C_col, ratio;
    int *current_D_block_rows, *local_ring_CprodD, *elements_CprodD, **CprodD;
    int *local_C, *local_D, *local_CplusD, *local_CprodB, *local_A;
    MPI_Request req[2];
    MPI_Status status[2];

    // Check if MPI environment got initialized correctly
    rc = MPI_Init(&argc, &argv);
    setbuf(stdout, NULL);
    if (rc != 0) {  // Else terminate MPI environment
        printf("MPI initialization failed...\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &p);  // Calculates number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // Gives rank to every process

    while(1) {
        // Master process 0 prints a menu of choices for user
        if (rank == 0) {
            printf("\n\nMENU\n\n1. Continue\n2. Exit\n\nChoice: ");
            fflush(stdout);
            scanf("%d", &choice);
        }
            // Master sends to every other process user' s choice
        MPI_Bcast(&choice, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (choice == 2) {
            if (rank == 0) {
                printf("\n\nProgramm terminated...\n");
            }

            MPI_Finalize();
            return 0;
        }


        if (rank == 0) {
            // Make sure that N%p = 0
            do {
                printf("\nType size N of arrays (A(1xN), B(Nx1), C(NxN) and D(NxN)): ");
                scanf("%d", &N);
            } while (N % p != 0);

            n = N/p;    // The number of rows every process will take on to calculate the final result
        }

        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);


        if (rank == 0) {
            // Allocate memory for each array
            A = (int*)malloc(N * sizeof(int));
            if (!A) {
                printf("Memory allocation error for array A...");
                MPI_Abort(MPI_COMM_WORLD, 100);
            }

            B = (int*)malloc(N * sizeof(int));
            if (!B) {
                printf("Memory allocation error for array B...");
                MPI_Abort(MPI_COMM_WORLD, 100);
            }

            C = (int**)malloc(N * sizeof(int*));
            // elements_C pointer is used in order for the allocated memory of EVERY row
            // of array C to be continuous in memory
            elements_C = (int*)malloc(N * N * sizeof(int));
            if (!C || !elements_C) {
                printf("Memory allocation error for array C...");
                MPI_Abort(MPI_COMM_WORLD, 100);
            }
            for (i = 0; i < N; i++) {
                C[i] = &(elements_C[i * N]);
            }

            D = (int**)malloc(N * sizeof(int*));
            elements_D = (int*)malloc(N * N * sizeof(int));
            if (!D || !elements_D) {
                printf("Memory allocation error for array D...");
                MPI_Abort(MPI_COMM_WORLD, 100);
            }
            for (i = 0; i < N; i++) {
                D[i] = &(elements_D[i * N]);
            }


            // Array CplusD is a N*N array and the result array of C + D
            CplusD = (int**)malloc(N * sizeof(int*));
            elements_CplusD = (int*)malloc(N * N * sizeof(int));
            if (!CplusD || !elements_CplusD) {
                printf("Memory allocation error for array CplusD...");
                MPI_Abort(MPI_COMM_WORLD, 100);
            }
            for (i = 0; i < N; i++) {
                CplusD[i] = &(elements_CplusD[i * N]);
            }

            // Array CprodB is a N*1 array and the result of C*B
            CprodB = (int*)malloc(N * sizeof(int));
            if (!CprodB) {
                printf("Memory allocation error for array CprodB...");
                MPI_Abort(MPI_COMM_WORLD, 100);
            }

            // Array CprodD is a N*N array and the result array of C * D
            CprodD = (int**)malloc(N * sizeof(int*));
            elements_CprodD = (int*)malloc(N * N * sizeof(int));
            if (!CprodD || !elements_CprodD) {
                printf("Memory allocation error for array CprodD...");
                MPI_Abort(MPI_COMM_WORLD, 100);
            }
            for (i = 0; i < N; i++) {
                CprodD[i] = &(elements_CprodD[i * N]);
            }


            // Get data for array A
            for (i = 0; i < N; i++) {
                printf("\nType element A1%d: ", i+1);
                scanf("%d", &A[i]);
            }

            // Get data for array B
            for (i = 0; i < N; i++) {
                printf("\nType element B%d1: ", i+1);
                scanf("%d", &B[i]);
            }

            // Get data for array C
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    printf("\nType element C%d%d: ", i+1, j+1);
                    scanf("%d", &C[i][j]);
                }
            }

            // Get data for array D
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    printf("\nType element D%d%d: ", i+1, j+1);
                    scanf("%d", &D[i][j]);
                }
            }


            // Prints every array clearly
            // Print array A
            printf("Array A:\n");
            for (i = 0; i < N; i++) {
                printf("%-4d  ", A[i]);
            }
            printf("\n");

            // Print array B
            printf("Array B:\n");
            for (i = 0; i < N; i++) {
                printf("%-4d\n", B[i]);
            }
            printf("\n");


            // Print array C
            printf("Array C:\n");
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    printf("%-4d  ", C[i][j]);
                }
                printf("\n");
            }
            printf("\n");

            // Print array D
            printf("Array D:\n");
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    printf("%-4d  ", D[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }


        //------------------------------------------------------------------------------------


        // I. C(NxN) + D(NxN)
        // Allocate memory for the sub-vector of each process
        local_C = (int*)malloc(n * N * sizeof(int));
        local_D = (int*)malloc(n * N * sizeof(int));
        local_CplusD = (int*)malloc(n * N * sizeof(int));
        if (!local_C || !local_D || !local_CplusD) {
            printf("Memory allocation error for local sub-vectors...");
            MPI_Abort(MPI_COMM_WORLD, 100); // Terminate EVERY process and return code 100
        }

        // Every process, including root (master),
        MPI_Scatter(elements_C, n*N, MPI_INT, local_C, n*N, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(elements_D, n*N, MPI_INT, local_D, n*N, MPI_INT, 0, MPI_COMM_WORLD);

        //Calculate result of local_C + local_D
        for (i = 0; i < n*N; i++) {
            local_CplusD[i] = local_C[i] + local_D[i];
        }

        // After every process, including root, calculates local_CplusD they are going to send their results to root (master)
        MPI_Gather(local_CplusD, n * N, MPI_INT, elements_CplusD, n * N, MPI_INT, 0, MPI_COMM_WORLD);


        //------------------------------------------------------------------------------------


        // II. C(NxN) * B(Nx1)
        // Array C is already sent to every process (local_C) so only array B remains to be sent
        if (rank != 0) {
            B = (int*)malloc(N * sizeof(int));
            if (!B) {
                printf("Memory allocation error for array B...");
                MPI_Abort(MPI_COMM_WORLD, 100);
            }
        }

        // Allocate memory for local_CprodB
        local_CprodB = (int*)malloc(n * sizeof(int));
        if (!local_CprodB) {
            printf("Memory allocation error for array local_CprodB...");
            MPI_Abort(MPI_COMM_WORLD, 100);
        }

        MPI_Bcast(B, N, MPI_INT, 0, MPI_COMM_WORLD);

        // Calculation: Every row of C gets multipled by vector B resulting to a N*1 array
        for (i = 0; i < n; i++) {
            local_CprodB[i] = 0;
            for (j = 0; j < N; j++) {
                local_CprodB[i] += local_C[i*N+j] * B[j];
            }
        }

        // Root collects the results of all processes including itself and saves them to CprodB
        MPI_Gather(local_CprodB, n, MPI_INT, CprodB, n, MPI_INT, 0, MPI_COMM_WORLD);

        free(local_CprodB);


        //------------------------------------------------------------------------------------


        // III. A(1xN) * B(Nx1)
        // Now array B is already sent to every process so array A will be sent devided
        local_A = (int*)malloc(n * sizeof(int));
        if (!local_A) {
            printf("Memory allocation error for array local_A...");
            MPI_Abort(MPI_COMM_WORLD, 100);
        }

        // Every process receives a part of A and saves it to local_A
        MPI_Scatter(A, n, MPI_INT, local_A, n, MPI_INT, 0, MPI_COMM_WORLD);

        local_AprodB = 0;
        for (i = 0; i < n; i++) {
            local_AprodB += local_A[i] * B[(rank*n)+i]; // Because B is whole, it must be devided manually
        }

        AprodB = 0;

        // Every process sends its local_AprodB to root (0) which is summing all of their results to 1 final result
        MPI_Reduce(&local_AprodB, &AprodB, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        free(local_A);


        //------------------------------------------------------------------------------------


        // IV. C(NxN) * D(NxN)

        // This block of code will be executed ONLY if N = p
        // Allocate memory for the local results of every process
        local_ring_CprodD = (int*)calloc(n * N, sizeof(int));   // calloc is used in order to initialize values of array local_ring_CprodD with 0
        // current_D_block_rows is a block of rows of D that is transfered between processes
        current_D_block_rows = (int*)malloc(n * N * sizeof(int));
        if (!local_ring_CprodD || !current_D_block_rows) {
            printf("Memory allocation error for array local_ring_CprodD or current_D_block_rows...");
            MPI_Abort(MPI_COMM_WORLD, 100);
        }

        source_p = (rank + 1) % p;  // Source process: It will send D data to its PREVIOUS process (dest process)
        dest_p = (rank - 1 + p) % p;    // Destination process: It will receive D data from source process

        localD_index = rank;

        for (i = 0; i < p; i++) {   // First for loop runs p times, as much as the number of processes, in order for every process to receive its local D data
            for (j = 0; j < n; j++) {   // Runs for every row of local_C that EVERY process owns
                for (k = 0; k < n; k++) {   // Runs for every row of local_D that the source process sent
                    C_col = (localD_index * n) + k; // C_col is column of C that must be multiplied by local_D

                    ratio = local_C[(j*N) + C_col]; // ratio is factor from local_C

                    // Sum multiplied row of ratio, of local_D to the result
                    for (l = 0; l < N; l++) {
                        local_ring_CprodD[(j*N) + l] += ratio * local_D[(k*N) + l];
                    }
                }
            }

            // If the program is not in the last step
            if (i < p - 1) {
                // Receive current_D_block_rows from the source process and save it
                MPI_Irecv(current_D_block_rows, n*N, MPI_INT, source_p, 0, MPI_COMM_WORLD, &req[0]);

                // Send local_D to destination process (previous process)
                MPI_Isend(local_D, n*N, MPI_INT, dest_p, 0, MPI_COMM_WORLD, &req[1]);

                // Wait unil the exchange has been completed
                MPI_Waitall(2, req, status);

                // Switch pointers
                temp = local_D;
                local_D = current_D_block_rows;
                current_D_block_rows = temp;

                // Update index
                localD_index = (localD_index + 1) % p;
            }


        }

        // Root collects the results of all processes including itself and saves them to CprodD
        MPI_Gather(local_ring_CprodD, n*N, MPI_INT, elements_CprodD, n*N, MPI_INT, 0, MPI_COMM_WORLD);


        //------------------------------------------------------------------------------------


        // Free allocated memory for sub-vectors
        free(local_ring_CprodD);
        free(current_D_block_rows);
        free(local_C);
        free(local_D);
        free(local_CplusD);
        free(B);

        if (rank == 0) {
            //I
            // Print result array CplusD of operation C + D
            printf("\n---------------------\nArray C + D:\n");
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    printf("%-4d  ", CplusD[i][j]);
                }
                printf("\n");
            }



            //II
            // Print result array CprodB of operation C * B
            printf("\n---------------------\nArray C * B:\n");
            for (i = 0; i < N; i++) {
                printf("%-4d\n", CprodB[i]);
            }



            // III
            // Print result, which is a number, of operatrion A * B
            printf("\n---------------------\n");
            printf("A * B: %d\n", AprodB);



            // IV
            // Print result, which is a N*N array, of operatrion C * D
            printf("\n---------------------\nArray C * D:\n");
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    printf("%-4d  ", CprodD[i][j]);
                }
                printf("\n");
            }
            printf("---------------------\n");


            // Free allocated memory for arrays
            free(A);

            free(elements_C);
            free(C);

            free(elements_D);
            free(D);

            free(elements_CplusD);
            free(CplusD);

            free(CprodB);
        }


        // In order to synchronize the outputs of every process
        MPI_Barrier(MPI_COMM_WORLD);

    }

    MPI_Finalize(); // Terminate MPI environment

    return 0;
}
