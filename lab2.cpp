#include <stdio.h>
#include <string.h>
#include <random>
#include <iostream>
#include <random>
#include <cmath>
#include <mpi.h>
using namespace std;

// Function to calculate the value of the integrand
double integrand(double x) {
    return sqrt(1 - x * x);
}

// Trapezoidal rule for numerical integration
double trapezoidal_rule(int n) {
    double a = 0.0;  // Lower limit of integration
    double b = 1.0;  // Upper limit of integration
    double h = (b - a) / n;
    double sum = 0.5 * (integrand(a) + integrand(b));

    for (int i = 1; i < n; ++i) {
        double x_i = a + i * h;
        sum += integrand(x_i);
    }

    return sum * h;
}

int main(int argc, char* argv[]){
   
   int N = 200;
   int x1 = 0;
   int x2 = 1;

   int subinterval = 10000;  // Number of subintervals (increase for higher accuracy)

   //This is for mpi
   int rank, size;
   double message [100];
   MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2) {
        cout << "This program requires at least 2 processes." << std::endl;
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        // Process 0 sends messages to all other processes
        for (int dest = 1; dest < size; dest++) {
            random_device rd;// Initialize with a random seed from the hardware
            mt19937 gen(rd());// Mersenne Twister pseudo-random number generator
            uniform_real_distribution<double> distribution(0.0, 1.0);// Create a distribution that generates numbers between 0 and 1
            double random_number = distribution(gen);// Generate a random number between 0 and 1
            MPI_Send(&random_number, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
        }

        // Process 0 receives messages from all other processes
        double total_sum = 0;
        for (int source = 1; source < size; source++) {
            double received_message;
            MPI_Recv(&received_message, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum += received_message;
            
        }
        double area = (total_sum / (size - 1)) * (x2 - x1);
        cout << "Empirical result: " << area << endl;
        double result = trapezoidal_rule(subinterval);
        cout << "Analytic Result: " << result << endl;

        } 
        
        else {
        // All other processes receive messages from Process 0
        double received_message;
        MPI_Recv(&received_message, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Send a message back to Process 0
        double sum = 0;
        sum = sqrt(1 - (received_message * received_message));
        MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }  
   
    MPI_Finalize();
    return 0;
}