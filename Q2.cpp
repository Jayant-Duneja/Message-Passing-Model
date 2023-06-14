#include <iostream>
#include <mpi.h>
#include <fstream>
#include <bits/stdc++.h>
using namespace std;

int main(int argc, char **argv)
{
    int rank, numprocs;

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
    vector<vector<int>> edges;
    vector<int> allvals;
    if (rank == 0)
    {

        ifstream inpfile(argv[1]);
        ofstream outfile(argv[2]);
        int v, e;
        int word_no = -2;
        while (inpfile >> txt)
        {
            // cout << txt << endl;
            if (word_no == -2)
            {
                v = stoi(txt);
                for (int i = 0; i < v; i++)
                {
                    vector<int> node;
                    for (int j = 0; j < v; j++)
                    {
                        node.push_back(-1);
                    }
                    edges.push_back(node);
                }
            }
            else if (word_no == -1)
            {
                e = stoi(txt);
            }
            else
            {
                int e1, e2, w;
                if (word_no % 3 == 0)
                {
                    e1 = stoi(txt);
                }
                if (word_no % 3 == 1)
                {
                    e2 = stoi(txt);
                }
                if (word_no % 3 == 2)
                {
                    w = stoi(txt);
                    edges[e1 - 1][e2 - 1] = w;
                    edges[e2 - 1][e1 - 1] = w;
                }
            }
            word_no++;

            // Output the text from the file
            //     num = stoi(txt);
        }
        MPI_Bcast(&v, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // cout << "Reached 1" << endl;
        int *adj = (int *)malloc((v * v) * sizeof(int));
        int k = 0;
        for (int i = 0; i < v; i++)
        {
            for (int j = 0; j < v; j++)
            {
                adj[k] = edges[i][j];
                k++;
            }
        }

        // for(int i = 0; i < v*v; i+= v){
        //             for(int j = 0; j < v; j++){
        //     cout << adj[i + j] << " ";
        // }
        // cout << endl;

        // }

        // All pairs of 3 vertices
        for (int i = 1; i <= v; i++)
        {
            for (int j = i + 1; j <= v; j++)
            {
                for (int k = j + 1; k <= v; k++)
                {
                    allvals.push_back(-1);
                    allvals.push_back(i);
                    allvals.push_back(j);
                    allvals.push_back(k);
                }
            }
        }
        // All pairs of 4 vertices
        for (int i = 1; i <= v; i++)
        {
            for (int j = i + 1; j <= v; j++)
            {
                for (int k = j + 1; k <= v; k++)
                {
                    for (int l = k + 1; l <= v; l++)
                    {
                        // cout << l << endl;
                        allvals.push_back(i);
                        allvals.push_back(j);
                        allvals.push_back(k);
                        allvals.push_back(l);
                    }
                }
            }
        }

        // cout << "Reached 2" << endl;

        // cout << "Reached2.1" << endl;
        int comb = allvals.size();
        // cout << "Reached2.1" << endl;
        // cout << comb << endl;
        int *sendvals = (int *)malloc((comb) * sizeof(int));
        // cout << "Reached 2.12" << endl;

        for (int i = 0; i < comb; i++)
        {
            sendvals[i] = allvals[i];
        }
        // cout << "Reached 2.13" << endl;

        int n = comb / (4 * numprocs);
        // cout << "Reached 2.2" << endl;

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(adj, v * v, MPI_INT, 0, MPI_COMM_WORLD);
        // cout << "Reached 2.3" << endl;
        int counts_send[numprocs] = {0}, displacements[numprocs] = {0};
        counts_send[0] += ((comb / 4) % numprocs) * 4;
        for (int i = 0; i < numprocs; i++)
        {
            counts_send[i] += n * 4;
            if (i > 0)
            {
                displacements[i] = displacements[i - 1] + counts_send[i - 1];
            }
        }
        // cout << "Reached 2.4" << endl;
        int *check = (int *)malloc((counts_send[0]) * sizeof(int));
        // cout << "Reached 3" << endl;

        // cout << "Scattering dims:" << endl;
        // for(int i = 0; i < numprocs; i++){
        //     cout << counts_send[i] << " ";
        // }
        // cout << endl;
        // for(int i = 0; i < numprocs; i++){
        //     cout << displacements[i] << " ";
        // }
        // cout << endl;
        MPI_Scatterv(sendvals, counts_send, displacements, MPI_INT, check, counts_send[0], MPI_INT, 0, MPI_COMM_WORLD);

        int l = 0;
        int r = counts_send[0];
        int dist[11] = {0};
        int x = 0;
        // cout << "Reached 4" << endl;

        for (int i = l; i < r; i += 4)
        {
            int edge_count = 0;
            if (check[i] == -1)
            {
                int temp[] = {check[i + 1], check[i + 2], check[i + 3]};
                // cout << x << " " << temp[0] << " " << temp[1] << " " << temp[2] << endl;
                int flag = 0;
                for (int j = 0; j < 3; j++)
                {
                    for (int k = j + 1; k < 3; k++)
                    {
                        if (adj[v * (temp[j] - 1) + temp[k] - 1] == -1)
                        {
                            flag = 1;
                            edge_count = -1;
                            break;
                        }
                        edge_count += adj[v * (temp[j] - 1) + temp[k] - 1];
                    }
                    if (flag == 1)
                    {
                        break;
                    }
                }
                if (flag == 0 && edge_count != -1)
                {
                    dist[edge_count] += 1;
                }
                // cout << edge_count << endl;
            }
            else
            {
                int temp[] = {check[i], check[i + 1], check[i + 2], check[i + 3]};
                // cout << x << " " << temp[0] << " " << temp[1] << " " << temp[2]  << " " << temp[3] << endl;

                int flag = 0;
                for (int j = 0; j < 4; j++)
                {
                    for (int k = j + 1; k < 4; k++)
                    {
                        if (adj[v * (temp[j] - 1) + temp[k] - 1] == -1)
                        {
                            flag = 1;
                            edge_count = -1;
                            break;
                        }
                        edge_count += adj[v * (temp[j] - 1) + temp[k] - 1];
                    }
                    if (flag == 1)
                    {
                        break;
                    }
                }
                if (flag == 0 && edge_count != -1)
                {
                    dist[4 + edge_count] += 1;
                }
                // cout << edge_count << endl;
            }
            x++;
        }

        // cout << "Reached 5" << endl;

        int alld[11 * numprocs];
        displacements[0] = 0;
        for (int i = 0; i < numprocs; i++)
        {
            counts_send[i] = 11;
            if (i > 0)
            {
                displacements[i] = displacements[i - 1] + counts_send[i - 1];
            }
        }
        MPI_Gather(dist, 11, MPI_INT, alld, 11, MPI_INT, 0, MPI_COMM_WORLD);

        for (int i = 0; i < 11; i++)
        {
            int total = 0;
            for (int j = 0; j < numprocs; j++)
            {
                total += alld[11 * j + i];
            }
            if (i < 4)
            {
                outfile << "3 " << i << " " << total << endl;
            }
            else
            {
                outfile << "4 " << i - 4 << " " << total << endl;
            }
        }
        outfile.close();
    }
    else
    {
        int v, n;
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // cout << n  << " rank " << rank << endl;

        MPI_Bcast(&v, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // cout << v  << " rank " << rank << endl;
        int *adj = (int *)malloc((n * n) * sizeof(int));
        // cout << v  << " rank " << rank << endl;
        MPI_Bcast(adj, n * n, MPI_INT, 0, MPI_COMM_WORLD);
        // cout << n  << " rank " << rank << endl;
        int *check = (int *)malloc((v * 4) * sizeof(int));
        MPI_Scatterv(NULL, NULL, NULL, MPI_INT, check, v * 4, MPI_INT, 0, MPI_COMM_WORLD);
        // cout << "received " << v << " by rank " << rank << endl;
        // for(int i = 0; i < v*4; i++)
        //     cout << check[i] << " ";
        // cout << endl;

        int l = 0;
        int r = v * 4;
        int dist[11] = {0};
        // cout << "Reached 6 " << r << endl;

        for (int i = l; i < r; i += 4)
        {
            // cout << check[i] << " " <<  check[i + 1] << " " <<  check[i + 2] << " " <<  check[i + 3]  << endl;
            int edge_count = 0;
            if (check[i] == -1)
            {
                int temp[] = {check[i + 1], check[i + 2], check[i + 3]};
                int flag = 0;
                for (int j = 0; j < 3; j++)
                {
                    for (int k = j + 1; k < 3; k++)
                    {
                        if (adj[n * (temp[j] - 1) + temp[k] - 1] == -1)
                        {
                            flag = 1;
                            edge_count = -1;
                            break;
                        }
                        edge_count += adj[n * (temp[j] - 1) + temp[k] - 1];
                    }
                    if (flag == 1)
                    {
                        break;
                    }
                }
                if (flag == 0 && edge_count != -1)
                {
                    dist[edge_count] += 1;
                }
            }
            else
            {
                int temp[] = {check[i], check[i + 1], check[i + 2], check[i + 3]};
                int flag = 0;
                for (int j = 0; j < 4; j++)
                {
                    for (int k = j + 1; k < 4; k++)
                    {
                        if (adj[n * (temp[j] - 1) + temp[k] - 1] == -1)
                        {
                            flag = 1;
                            edge_count = -1;
                            break;
                        }
                        edge_count += adj[n * (temp[j] - 1) + temp[k] - 1];
                    }
                    if (flag == 1)
                    {
                        break;
                    }
                }
                if (flag == 0 && edge_count != -1)
                {
                    dist[4 + edge_count] += 1;
                }
            }
            // cout << "exit time" << endl;
        }
        // cout << "Reached 7" << endl;

        MPI_Gather(dist, 11, MPI_INT, NULL, 0, NULL, 0, MPI_COMM_WORLD);
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