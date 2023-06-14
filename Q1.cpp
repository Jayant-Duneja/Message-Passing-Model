#include <iostream>
#include <mpi.h>
#include <fstream>
#include <bits/stdc++.h>
using namespace std;

int main(int argc, char **argv)
{
    int rank, numprocs;
    ifstream inpfile(argv[1]);
    ofstream outfile(argv[2]);
    // initiate MPI
    MPI_Init(&argc, &argv);

    // get size of the current communicator
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    // get current process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* synchronize all processes */
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    // enter your code here
    
    string txt;
    int num;
    if (rank == 0)
    {
        while (getline(inpfile, txt))
        {
            // Output the text from the file
            num = stoi(txt);
        }

        int prime = 1;
        int crank = rank + 2;
        while (crank <= sqrt(num) && prime == 1)
        {
            if (num % crank == 0)
            {
                prime = 0;
                break;
            }
            crank += numprocs;
        }
        int is_prime;
        for (int i = 1; i < numprocs; i++)
        {
            MPI_Send(&num, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Recv(&is_prime, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (prime == 0 || is_prime == 0)
            {
                prime = 0;
            }
        }

        if (prime == 1)
        {
            outfile << "YES";
        }
        else
        {
            outfile << "NO";
        }
        outfile.close();
    }
    else
    {
        int num;
        int prime = 1;
        int crank = rank + 2;
        MPI_Recv(&num, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (num <= 3)
        {
            MPI_Send(&prime, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        else
        {
            while (crank <= sqrt(num) && prime == 1)
            {
                if (num % crank == 0)
                {
                    prime = 0;
                }
                crank += numprocs;
            }
            MPI_Send(&prime, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime() - start_time;
    double maxTime;
    // get max program run time for all processes
    MPI_Reduce(&end_time, &maxTime, 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        cout << "Total time (s): " << maxTime << "\n";
    }
    // shut down MPI and close
    MPI_Finalize();
    return 0;
}