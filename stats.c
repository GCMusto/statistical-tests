#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

// Define the structure for a vector
typedef struct {
  double* data;
  size_t size;
} Vector;

// Define structure for holding the test
typedef struct {
  Vector* data1;              // pointer to sample 1
  Vector* data2; 	      // pointer to sample 2
  size_t size1;		      // size of the first sample
  size_t size2;		      // size of the second sample
  Vector* ranks1;	      // pointer to vector of ranks for sample 1
  Vector* ranks2;	      // pointer to vector of ranks for sample 2
  double u1;		      // Mann-Whitney U statistics for sample 1
  double u2;		      // Mann-Whitney U statistics for sample 2
  double U;                   // Mann-Whitney U statistics
  double mu;                  // Mann-Whitney U statistics mean
  double sd;                  // Mann-Whitney U statistics standard deviation
  const char* alternative;    // two-sided, less, greater
  bool correction;            // true for continuity correction
  double pvalue;	      // p-value corresponding to Mann-Whitney test
} MannWhitney_data;

// Function to determine the sign of a double
int sign(double num) {
    if (num > 0) {
        return 1; // Positive
    } else if (num < 0) {
        return -1; // Negative
    } else {
        return 0; // Zero
    }
}

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

void split_ranking(Vector* ranked_vector, Vector* a, Vector* b) {

 for (size_t i = 0; i < a->size; i++) {
    a->data[i] = ranked_vector->data[i];
  }

  for (size_t i = 0; i < b->size; i++) {
    b->data[i] = ranked_vector->data[a->size + i];
  }

}


double mwmu(size_t n1, size_t n2)  {
  return ((double)n1 * (double)n2)/2.0;
}

double mwsd(Vector* rank1, Vector* rank2) {

    double tie_adjustment = 0.0;
    size_t n1 = rank1->size;
    size_t n2 = rank2->size;    
    size_t N = n1 + n2;
    Vector rank = combine_vectors(rank1, rank2);
    
    for (size_t i = 0; i < N; i++) {
        double current_rank = rank.data[i];
        size_t tie_count = 0;

        // Count the number of tied values for the current rank
        for (size_t j = 0; j < N; j++) {
            if (rank.data[j] == current_rank) {
                tie_count++;
            }
        }

        // Calculate the difference between the cube of ties and ties for distinct ranks
        tie_adjustment += (tie_count * tie_count * tie_count) - tie_count;
    }

    double sd = sqrt(((n1 * n2)/12.0) * ((N + 1) - (tie_adjustment/(N * (N - 1)))));
    return sd;
  
}

// Function to calculate the p-value for a two-tailed normal distribution
double calculate_pvalue(double observedValue, double mean,
			double stdDev, const char* alternative,
			bool correction) {

    double z = 0.0;
    double p_value = 0.0;
			
    if (strcmp(alternative, "two-sided") == 0) {
      z = (observedValue - mean - correction * sign(observedValue - mean) * 0.5)/stdDev;      
      p_value = 2.0 * (1.0 - 0.5 * (1.0 + erf(z / sqrt(2.0))));
    } else if (strcmp(alternative, "less") == 0) {
      z = (observedValue - mean + correction * 0.5)/stdDev;
      p_value = 1.0 - 0.5 * (1.0 + erf(z / sqrt(2.0)));
    } else if (strcmp(alternative, "greater") == 0) {
      z = (observedValue - mean - correction * 0.5)/stdDev;
      p_value = 0.5 * (1.0 + erf(z / sqrt(2.0)));
    } else {
      fprintf(stderr, "Wrong alternative hypothesis selected!\n");
      exit(1);
    }
    
//    if (z < 0) {
//        // Calculate the left tail of the distribution
//        p_value = 0.5 * (1.0 + erf(z / sqrt(2.0)));
//    } else {
//        // Calculate the right tail of the distribution
//        p_value = 0.5 * (1.0 - erf(-z / sqrt(2.0)));
//    }
//
    return p_value;
}

void mwtest(MannWhitney_data* mwData){

  mwData->u1 = 0;
  mwData->u2 = 0;  
  size_t n1 = mwData->size1;
  size_t n2 = mwData->size2;
  double sum_rank1 = 0;
  double sum_rank2 = 0;

  
  for (size_t i = 0; i < n1; i++){
    sum_rank1 += mwData->ranks1->data[i];
  }

  sum_rank2 = (n1 + n2) * ((n1 + n2) + 1.0)/2.0;
  //mwData->u1 = (n1 * n2) + (n1 * (n1 + 1)/2) - sum_rank1;
  //mwData->u2 = (n1 * n2) + (n2 * (n2 + 1)/2) - sum_rank2;
  mwData->u1 = sum_rank1 - (n1 * (n1 + 1))/2.0;
  mwData->u2 = sum_rank2 - (n2 * (n2 + 1))/2.0;

  if (mwData->u1 > mwData->u2) mwData->U = mwData->u2;
  else mwData->U = mwData->u1;

  mwData->pvalue =
    calculate_pvalue(
      mwData->U,
      mwData->mu,
      mwData->sd,
      mwData->alternative,
      mwData->correction);
  
}

void print_mw_test(MannWhitney_data* mwData) {

  printf("========================================================================\n");
  printf("=                                                                      \n");  
  printf("= data: vector1 and vector2                                            \n");
  printf("= U = %.2f, p-value = %.4f                                             \n", mwData->U, mwData->pvalue);
  printf("= Alternative hypothesis: true location shift is not equal to zero     \n");
  printf("= \n");
  printf("= U1 = %.2f, U2 = %.2f, mu = %.2f, sd = %.4f                           \n", mwData->u1, mwData->u2, mwData->mu, mwData->sd);
  printf("=                                                                      \n");
  printf("========================================================================\n");
  
}

int main(){
  size_t size1, size2;
    
  // Read the necessary data
  printf("== Enter the size of vector 1: \n");
  scanf("%zu", &size1);
  printf("\n== Enter the size of vector 2: \n");
  scanf("%zu", &size2);
  
  Vector vector1 = create_vector(size1);
  Vector vector2 = create_vector(size2);
  
  printf("\n\n== Vector 1 ========================\n");
  input_vector(&vector1);
  printf("\n== Vector 2 ========================\n");  
  input_vector(&vector2);
  printf("\n\n");
  
  // Perform observation ranking
  Vector ab = combine_vectors(&vector1, &vector2);
  Vector ranked_vector = rank_vector(&ab, true);

  // Split ranking in two vectors
  Vector rank1 = create_vector(size1);
  Vector rank2 = create_vector(size2);
  split_ranking(&ranked_vector, &rank1, &rank2);

  // Populate Mann-Whitney test structure
  MannWhitney_data mwData;
  mwData.data1 = &vector1;
  mwData.data2 = &vector2;
  mwData.size1 = vector1.size;
  mwData.size2 = vector2.size;
  mwData.ranks1 = &rank1;
  mwData.ranks2 = &rank2;
  mwData.mu = mwmu(mwData.size1, mwData.size2);
  mwData.sd = mwsd(mwData.ranks1, mwData.ranks2);  
  mwData.alternative = "two-sided";
  mwData.correction = true;
  mwtest(&mwData);
  print_mw_test(&mwData);
  
  // Free resources
  free_vector(&vector1);
  free_vector(&vector2);
  free_vector(&ab);
  free_vector(&ranked_vector);
  free_vector(&rank1);
  free_vector(&rank2);
  
  return 0;
}
