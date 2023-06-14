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
    int n;
    if (rank == 0)
    {
        inpfile >> txt;
        n = stoi(txt);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double *matrix = (double *)malloc((n * (n * 1)) * sizeof(double));
    double ans[n];

    if (rank == 0)
    {
        for (int i = 0; i < n * (n + 1); i++)
        {
            inpfile >> matrix[i];
        }
    }
    MPI_Bcast(matrix, n * (n + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < n - 1; i++)
    {
        // cout << "Entry for rank " << rank << endl;
        // for (int j = 0; j < n * (n + 1); j++)
        // {
        //     cout << matrix[j] << " ";
        // }
        // cout << endl;
        if (rank == 0)
        {
            double cur_max = abs(matrix[(n + 1) * i + i]);
            int max_ind = i;
            for (int j = i + 1; j < n; j++)
            {
                if (abs(matrix[j * (n + 1) + i]) > cur_max)
                {
                    cur_max = abs(matrix[j * (n + 1) + i]);
                    max_ind = j;
                }
            }
            if (max_ind != i)
            {
                for (int j = 0; j <= n; j++)
                {
                    double temp = matrix[i * (n + 1) + j];
                    matrix[i * (n + 1) + j] = matrix[max_ind * (n + 1) + j];
                    matrix[max_ind * (n + 1) + j] = temp;
                }
            }
        }
        // cout << "After Swap for ";
        // for (int j = 0; j < n * (n + 1); j++)
        // {
        //     cout << matrix[j] << " ";
        // }
        // cout << endl;
        int sending = n - i - 1;
        int counts_send[numprocs] = {0}, displacements[numprocs] = {0};
        counts_send[0] += ((sending) % numprocs) * (n + 1);
        for (int j = 0; j < numprocs; j++)
        {
            counts_send[j] += (sending / numprocs) * (n + 1);
            if (j > 0)
            {
                displacements[j] = displacements[j - 1] + counts_send[j - 1];
            }
        }
        double *cur_mat = (double *)malloc((counts_send[rank]) * sizeof(double));

        double cur_row[n + 1];
        if (rank == 0)
        {
            for (int j = 0; j < n + 1; j++)
            {
                cur_row[j] = matrix[(n + 1) * i + j];
            }
        }
        MPI_Bcast(cur_row, n + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {

            double *buf = (double *)malloc((sending * (n + 1)) * sizeof(double));
            for (int j = (i + 1) * (n + 1); j < n * (n + 1); j++)
            {
                buf[j - (i + 1) * (n + 1)] = matrix[j];
            }
            // cout << "Buff Enter: " << rank << "\n ";
            // for (int j = 0; j < sending * (n + 1); j++)
            // {
            //     cout << buf[j] << " ";
            // }
            // cout << endl;

            MPI_Scatterv(buf, counts_send, displacements, MPI_DOUBLE, cur_mat, counts_send[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            free(buf);
        }
        else
        {
            MPI_Scatterv(NULL, NULL, NULL, NULL, cur_mat, counts_send[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        int scale_ind = i;
        double scale_num = cur_row[scale_ind];
        // cout << "Cur Row init: " << rank << "\n ";
        // cout << cur_row[scale_ind]  << endl;
        for (int j = 0; j <= n; j++)
        {
            // cout << cur_row[j] << " " <<  cur_row[scale_ind] << " ";
            cur_row[j] = cur_row[j] / scale_num;
        }
        // cout << endl;
        // cout << "Cur Row: " << rank << "\n ";
        // for (int j = 0; j < (n + 1); j++)
        // {
        //     cout << cur_row[j] << " ";
        // }
        // cout << endl;
        double scale_val;
        for (int j = 0; j < counts_send[rank] / (n + 1); j++)
        {
            scale_val = cur_mat[(n + 1) * j + scale_ind];
            // cout << scale_val << " " << endl;
            for (int k = 0; k <= n; k++)
            {
                cur_mat[(n + 1) * j + k] -= scale_val * cur_row[k];
            }
        }
        // cout << "Cur Mat: " << rank << "\n ";
        // for (int j = 0; j < counts_send[rank]; j++)
        // {
        //     cout << cur_mat[j] << " ";
        // }
        // cout << endl;
        if (rank == 0)
        {
            double *buff = (double *)malloc((sending * (n + 1)) * sizeof(double));
            MPI_Gatherv(cur_mat, counts_send[rank], MPI_DOUBLE, buff, counts_send, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // cout << "Buff: " << rank << "\n ";
            // for (int j = 0; j < sending * (n + 1); j++)
            // {
            //     cout << buff[j] << " ";
            // }
            // cout << endl;
            for (int j = 0; j < sending * (n + 1); j++)
            {
                matrix[(i + 1) * (n + 1) + j] = buff[j];
            }
            // cout << "New Buf " << rank << "\n ";
            // for (int j = 0; j < n * (n + 1); j++)
            // {
            //     cout << matrix[j] << " ";
            // }
            // cout << endl;
            free(buff);
        }
        else
        {
            MPI_Gatherv(cur_mat, counts_send[rank], MPI_DOUBLE, NULL, NULL, NULL, NULL, 0, MPI_COMM_WORLD);
        }
        free(cur_mat);
    }

    if (rank == 0)
    {

        for (int i = n - 1; i >= 0; i--)
        {
            ans[i] = matrix[i * (n + 1) + n];
            for (int j = i + 1; j < n; j++)
            {
                ans[i] -= matrix[i * (n + 1) + j] * ans[j];
            }
            ans[i] /= matrix[i * (n + 1) + i];
        }

        for (int i = 0; i < n; i++)
            outfile << ans[i] << " ";
        outfile << endl;
    }

    outfile.close();

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
    free(matrix);

    // shut down MPI and close
    MPI_Finalize();

    return 0;
}