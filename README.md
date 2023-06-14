# Assignment 1
**Tanvi Kamble**\
**2018114004**

## Question 1:
- Considering the case of only a single process, we have to check if the number is divisble by any numbers in the range of 1 to square root of that number. 
- I have parallelised this block of code by making sure that each process only checks a subset of the numbers from 1 to Sqrt(n).
- Each rank starts from the number (2+rank) and checks for divisibility for all numbers upto Sqrt(n) with increments of the value of the number of processes.I then send the flag of each process and perform the logical AND operation with all the other flags from the other processes.
- If the final flag is true, the number is prime and otherwise the number is not prime.

## Question 2:
- In this question, first of all , I create the adjacency matrix for the given graph and broadcast that to all the processes.
- After this, I generated all possible combinations of size 3 and 4 of the nodes in the graph. I scattered these combinations to all the processes using the MPI_Scatterv function.
- Now, each process has a array of arrays where each element is a combination of the nodes of size 3 or 4. Each process checks if for each combination and increments the count of each type of clique in an array.
- After this, I gather all these count vectors using the MPI_Gatherv function and I display the counts after this.


## Question 3:

We partition the matrix and every process does the computations for their partition

- Process 0 reads input. Rows in the matrix are swapped to have the highest absolute for every column at the diagonal. This helps prevent computational errors due to floating point numbers.

- Every process iterates over every row
- If current row index lies in the range of processor then the row is normalized by dividing by pivot and broadcasted to the other process. The row is then subtracted by all rows after the current row index that belong to the current processor.
- If current row index lies above range of processor then the row is subtracted by all rows that belong to the current processor.
- This gives us the upper triangular form
- We now use back substitution to get the final answers
	- We iterate over the rows backwards.
	- Elements other than pivot element are eliminated by results of the previous rows
	- Answer is updated for pivot index and broadcasted to every process.
- Result is reported by process 0