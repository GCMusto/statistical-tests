#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// Define the structure for a vector
typedef struct {
  double* data;
  size_t size;
} Vector;

// Function to create a new vector
Vector create_vector(size_t size) {
  Vector vec;
  vec.size = size;
  vec.data = (double*)malloc(size * sizeof(double));
  if (vec.data == NULL) {
    fprintf(stderr, "Memory allocation failed!\n");
    exit(1);
  }
  return vec;
}

// Function to input values into a vector
void input_vector(Vector* vec){
  printf("Enter %zu values for the vector: \n", vec->size);
  for (size_t i = 0; i < vec->size; i++) {
    scanf("%lf", &vec->data[i]);
  }
}

// Free memory allocated for a vector
void free_vector(Vector* vec) {
  free(vec->data);
}

// Mann-Whitney test
// Function to combine two vectors
Vector combine_vectors(const Vector* a,const Vector* b) {

  Vector ab;
  ab.size = a->size + b->size;
  ab.data = (double*)malloc(ab.size * sizeof(double));
  if (ab.data == NULL) {
    fprintf(stderr, "Memory allocation failed!\n");
    exit(1);
  }

  // Copy elements from first vector
  for (size_t i = 0; i < a->size; i++) {
    ab.data[i] = a->data[i];
  }

  // Copy elements from first vector
  for (size_t i = 0; i < b->size; i++) {
    ab.data[a->size + i] = b->data[i];
  }
  
  return ab;
  
}

// Function to rank the values of a vector with handling for ties and assigning midpoint ranks
Vector rank_vector(Vector* vec, bool ascending) {

  Vector ranked_vector = create_vector(vec->size);
  		    		    
  for (size_t i = 0; i < vec->size; i++) {
    double value = vec->data[i];
    size_t tie_count = 1;
    size_t rank = 1;

    for (size_t j = 0; j < vec->size; j++) {
	if ((ascending && i != j && value > vec->data[j]) ||
	    (!ascending && i != j && value < vec->data[j])) {
	  rank++;
	}
	if ((ascending && i != j && value == vec->data[j]) ||
	    (!ascending && i != j && value == vec->data[j])) {
	  tie_count++;
	}
    }

    ranked_vector.data[i] = (double)rank + 0.5 * (double)(tie_count - 1);
    
  }

  return ranked_vector;
}


double mwmu(Vector* a, Vector* b)  {
  return (a->size * b->size)/2.0;
}

double mwsd(Vector* rank, Vector* a, Vector* b) {

    double tie_adjustment = 0.0;
    double n1n2 = a->size * b->size;
    double N = a->size * b->size;
    
    for (size_t i = 0; i < rank->size; i++) {
        double current_rank = rank->data[i];
        size_t tie_count = 0;

        // Count the number of tied values for the current rank
        for (size_t j = 0; j < rank->size; j++) {
            if (rank->data[j] == current_rank) {
                tie_count++;
            }
        }

        // Calculate the difference between the cube of ties and ties for distinct ranks
        tie_adjustment += (tie_count * tie_count * tie_count) - tie_count;
    }

    return sqrt((n1n2/12.0) * ((N + 1) - (tie_adjustment/(N * (N - 1)))));
  
}

// Function to calculate the p-value for a two-tailed normal distribution
double calculate_pvalue(double observedValue, double mean, double stdDev) {

    double z = (observedValue - mean) / stdDev;
    double p_value = 0.0;

    if (z < 0) {
        // Calculate the left tail of the distribution
        p_value = 0.5 * (1.0 + erf(z / sqrt(2.0)));
    } else {
        // Calculate the right tail of the distribution
        p_value = 0.5 * (1.0 - erf(-z / sqrt(2.0)));
    }

    return p_value;
}

double mwtest(Vector* rank, Vector* a, Vector* b){

  double sum_rank1 = 0;
  double sum_rank2 = 0;
  double u1 = 0;
  double u2 = 0;
  double mann_whitney_stat = 0;
  
  for (size_t i = 0; i < a->size; i++){
    sum_rank1 += rank->data[i];
  }

  sum_rank2 = (a->size + b->size) * ((a->size + b->size) + 1.0)/2.0;
  u1 = (a->size * b->size) + (a->size * (a->size + 1)/2) - sum_rank1;
  u2 = (a->size * b->size) + (b->size * (b->size + 1)/2) - sum_rank2;

  if (u1 > u2) mann_whitney_stat = u1;
  else mann_whitney_stat = u2;

  return calculate_pvalue(mann_whitney_stat, mwmu(a, b), mwsd(rank, a, b));
  
}

int main(){
  size_t size1, size2;

  // Read the necessary data
  printf("Enter the size of vector 1: \n");
  scanf("%zu", &size1);
  printf("Enter the size of vector 2: \n");
  scanf("%zu", &size2);

  Vector vector1 = create_vector(size1);
  Vector vector2 = create_vector(size2);

  printf("Vector 1 ========================\n");
  input_vector(&vector1);
  printf("Vector 2 ========================\n");  
  input_vector(&vector2);

  // Perform the test
  Vector ab = combine_vectors(&vector1, &vector2);
  Vector ordered_ab = rank_vector(&ab, true);
  
  for (size_t i = 0; i < ab.size; i++) {
    printf("%f \n", ordered_ab.data[i]);
  }  

  //  printf("\nThe Mann-Whitney statistic is: %f\n", mwtest(&ordered_ab, &vector1, &vector2));
  printf("\nThe Mann-Whitney p-value is: %f\n", mwtest(&ordered_ab, &vector1, &vector2));
  printf("\nThe Mann-Whitney statistic mean is: %f\n", mwmu(&vector1, &vector2));
  printf("\nThe Mann-Whitney statistic sd is: %f\n", mwsd(&ordered_ab, &vector1, &vector2));
    
  // Free resources
  free_vector(&vector1);
  free_vector(&vector2);
  free_vector(&ab);
  free_vector(&ordered_ab);
  
  return 0;
}
