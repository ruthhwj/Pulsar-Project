//START REGION RELEASE
//START REGION DEVELOP
/*
  The double precision functions are in fitting.c, and the other
  precision code (fitting_ld are long doubles) is derived from
  fitting.c. So always change the code in fitting.c and re-create the
  other precisions code. 

  For other precisions: Compile without math.h to see if all functions
  are in correct precision.

  double                              -> long double
  linalg_solve_matrix_eq_gauss_jordan -> linalg_solve_matrix_eq_gauss_jordan_ld
  fabs                                -> fabsl

 */


//START REGION RELEASE
#include <math.h>  // Take this line out to see what maths functions are used, i.e. all long double versions etc?
#include <stdlib.h>
#include <string.h>
#include "psrsalsa.h"
#include <gsl/gsl_linalg.h>
//START REGION DEVELOP


/*
  Given the n*n matrix a and collection of answer vectors b, find solutions vertors x of the matrix eq: a*x=b. Each vector b thus has n elements. In total there are m solution vectors b provided in a matrix, which is just a list of numbers. The m solution vectors are ordered as:

b[] = [b11...b1n,b21...b2n,...,bm1...bmn]

So the elements of each vector b appears consecutive in memory.

The matrix a is ordered such that the first n elements correspond to a row that is used to get the first element of the solution vectors.

The output matrix x has the same format as b (memory is never allocated and is only filled when findsolutions is nonzero.

The output matrix a_inv has the same format as a, and will contain the inverse of the matrix a. Memory should already be allocated when findinverse is nonzero, otherwise this matrix is not filled.

If preserve_a is zero, it is not guaranteed that the input array a is perserved (will not be the case when gsl is used).

If checkcaninverted is nonzero, a check is done if the matrix can be inverted. The check is always done in NR, but gsl returns by default nan's/inf's when it cannot be inverted.

Return values:
0 = ok
1 = memory allocation error
2 = Matrix a cannot be inverted
3 = Method not implemented
 */

int linalg_solve_matrix_eq(double *a, int n, int preserve_a, double *b, int m, double *x, int findsolutions, double *a_inv, int findinverse, int checkcaninverted, verbose_definition verbose)
{
  int usenr = 0;
  if(usenr) {
#ifdef NRAVAIL
    int i, j, k;
    long double **a_nr, **b_nr;
    a_nr = matrix_ld(1, n, 1, n);
    b_nr = matrix_ld(1, n, 1, m);
    if(a_nr == NULL || b_nr == NULL) {
      printerror(verbose.debug, "ERROR linalg_solve_matrix_eq: Cannot allocate memory.\n");
      return 1;
    }
    k = 0;
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++) {
	a_nr[1+i][1+j] = a[k++];
      }
    }
    k = 0;
    for(j = 0; j < m; j++) {
      for(i = 0; i < n; i++) {
	b_nr[1+i][1+j] = b[k++];
      }
    }
    
    if(gaussj_ld(a_nr, n, b_nr, m, 1) == 0) {  // Failed
      printwarning(verbose.debug, "WARNING gaussj failed: Matrix cannot be inverted.");
      free_matrix_ld(a_nr,1,n,1,n); 
      free_matrix_ld(b_nr,1,n,1,m); 
      return 2;
    }
    if(findinverse) {
      k = 0;
      for(i = 0; i < n; i++) {
	for(j = 0; j < n; j++) {
	  a_inv[k++] = a_nr[1+i][1+j];
	}
      }
    }
    if(findsolutions) {
      k = 0;
      for(j = 0; j < m; j++) {
	for(i = 0; i < n; i++) {
	  x[k++] = b_nr[1+i][1+j];
	}
      }
    }
    free_matrix_ld(a_nr,1,n,1,n); 
    free_matrix_ld(b_nr,1,n,1,m); 
#else
    printerror(verbose.debug, "ERROR linalg_solve_matrix_eq: No numerical recipies available when library was compiled.");
    return 3;
#endif
  }else {
    gsl_matrix_view a_matrix;
    gsl_permutation *p;
    double *a_copy;
    int j, sign;
    // The following function does not much, it only keeps the information about what the size of the array is. No memory is allocated and no memory is copied. So array a should not be freed until you're done with a_matrix.
    if(preserve_a) {
      a_copy = malloc(n*n*sizeof(double));
      if(a_copy == NULL) {
	printerror(verbose.debug, "ERROR linalg_solve_matrix_eq: Cannot allocate memory.\n");
	return 1;
      }
      memcpy(a_copy, a, n*n*sizeof(double));
    }
    a_matrix = gsl_matrix_view_array(a, n, n);
    p = gsl_permutation_alloc(n);
    // Find LU decomposition
    gsl_linalg_LU_decomp(&a_matrix.matrix, p, &sign);
    if(findinverse) {
      // Set inverse matrix to be at location of the output array
      gsl_matrix_view a_inv_matrix;
      a_inv_matrix = gsl_matrix_view_array(a_inv, n, n);
      // Fill the output array with the inverse of the matrix
      gsl_linalg_LU_invert(&a_matrix.matrix, p, &a_inv_matrix.matrix);
      if(checkcaninverted) {
	for(j = 0; j < n*n; j++) {
	  if(!isfinite(a_inv[j])) {
	    printwarning(verbose.debug, "WARNING gsl_linalg_LU_invert failed: Matrix cannot be inverted.");
	    gsl_permutation_free(p);
	    if(preserve_a) {
	      memcpy(a, a_copy, n*n*sizeof(double));
	      free(a_copy);
	    }
	    return 2;
	  }
	}
      }
    }
    if(findsolutions) {
      gsl_vector_view b_vector, x_vector;
      for(j = 0; j < m; j++) {
	// Select one vector b out of the input
	b_vector = gsl_vector_view_array(&b[j*n], n);
	// Select one vector x out of the output
	x_vector = gsl_vector_view_array(&x[j*n], n);
	// Obtain solution
	gsl_linalg_LU_solve(&a_matrix.matrix, p, &b_vector.vector, &x_vector.vector);
      }
      if(checkcaninverted) {
	for(j = 0; j < n*m; j++) {
	  if(!isfinite(x[j])) {
	    printwarning(verbose.debug, "WARNING gsl_linalg_LU_solve failed: Matrix cannot be inverted.");
	    gsl_permutation_free(p);
	    if(preserve_a) {
	      memcpy(a, a_copy, n*n*sizeof(double));
	      free(a_copy);
	    }
	    return 2;
	  }
	}
      }
    }
    gsl_permutation_free(p);
    if(preserve_a) {
      memcpy(a, a_copy, n*n*sizeof(double));
      free(a_copy);
    }
  }
  return 0;
}

//START REGION DEVELOP
//START REGION RELEASE

/* If returnerror is set, the singular matrix error which terminates
   program is captured and causes the function to return 0 rather than
   1. If not set, the program always returns 1, or terminates.

   Given the n*n matrix a and collection of answer vectors b, find solutions vertors x of the matrix eq: a*x=b. Each vector b thus has n elements. In total there are m solution vectors b provided in a matrix, which is just a list of numbers. The m solution vectors are ordered as:

   b[] = [b11...b1n,b21...b2n,...,bm1...bmn]

   So the elements of each vector b appears consecutive in memory.

   The matrix a is ordered such that the first n elements correspond to a row that is used to get the first element of the solution vectors.

   The input matrix matrixa (an n X n matrix) is replaced by its inverse.
   The input matrix matrixb (an n X m matrix) is replaced by the solution (same dimensions).

   Return 0 = Success
   Return 1 = Memory allocation error
   Return 2 = Singular matrix, no solution can be determined
 */
int linalg_solve_matrix_eq_gauss_jordan(double *matrixa, double *matrixb, int n, int m, verbose_definition verbose)
{
  int *orig_pivot_row, *orig_pivot_column, *pivot_column_done;
  int pivotnr, pivot_col, pivot_row, rowi, coli, i;
  double *ptr1, *ptr2, tmpvalue;
  double largest_element_value;
  int currow;

  orig_pivot_row = malloc(n*sizeof(int));
  orig_pivot_column = malloc(n*sizeof(int));
  pivot_column_done = malloc(n*sizeof(int));
  if(orig_pivot_row == NULL || orig_pivot_column == NULL || pivot_column_done == NULL) {
    printerror(verbose.debug, "ERROR linalg_solve_matrix_eq_gauss_jordan: Memory allocation error.");
    return 1;
  }
  for(coli = 0; coli < n; coli++)
    pivot_column_done[coli] = 0;    // Clear the array indicating which column has been processed

  // Loop over the total number of pivot points that need to be identified in matrix a???
  // Think of it as looping over the columns (in a non-sequential order), while the order of the rows will keep changing.
  for(pivotnr = 0; pivotnr < n; pivotnr++) {

    // Find a suitable "pivot element". The central idea of GJ
    // elimination is to imagine the matrix equation:
    // A*[x,Ainv] = [B,1]
    // Here [] means to paste matrices together and 1 is the identity matrix
    // x is the solution we want to find, which multiplied with A should give B.
    // The inverse of A, when multiplied by Ainv, should give the identity matrix.
    // The way to solve this is to realise that rows/columns can be changed 
    // in the matrix equation, or rows can be replaced by linear combinations 
    // of rows, without invalidating the equations (if all corresponding swaps 
    // are applied to all matrices). 
    // These operations can be applied to change A into the identity matrix, 
    // in which case the equation changed into
    // 1*[x,Ainv] = [x,Ainv]
    // Note that vectors x and the initial identity matrix do not explicity
    // exist in memory in the algorithm.
    largest_element_value = -1.0;  // Initialise max found value to be negativenegative, so the absolute values of the remaining elements always exceeds the initialisation value
    // Loop over all rows of the matrix
    for(rowi = 0; rowi < n; rowi++) {
      if(pivot_column_done[rowi] != 1) {  // Check if row wasn't selected=processed once before
	for(coli = 0; coli < n; coli++) {  // Loop over all columns of the not yet processed row
	  if(pivot_column_done[coli] == 0) {  // Check if rows=column? wasn't selected=processed before. I think if a given row is processed (1 on diagonal), then all values in the corresponding column will be set to zero. It doesn't quite make sense, as otherwise why does it need a check again? And surely the process needs to be looped over n times, implying each row need to be processed once, implying each pivot_column_done element is used one at a time.
	    tmpvalue = matrixa[rowi*n+coli];
	    if(fabs(tmpvalue) >= largest_element_value) {
	      largest_element_value = fabs(tmpvalue);
	      pivot_row=rowi;
	      pivot_col=coli;
	    }
	  }else if(pivot_column_done[coli] > 1) { //Not sure why this can happen. 
	    printerror(verbose.debug, "ERROR linalg_solve_matrix_eq_gauss_jordan: The matrix equation is singular, no solution can be determined.");
	    return 2;
	  }
	}
      }
    }
    // If the pivot element is zero, then the matrix 
    // cannot be re-written as the identity matrix. This means that
    // the solution does not exists.
    if(largest_element_value <= 0.0) {
      printerror(verbose.debug, "ERROR linalg_solve_matrix_eq_gauss_jordan: The matrix equation is singular, no solution can be determined.");
      return 2;
    }

    pivot_column_done[pivot_col] += 1; // Keep track of columns being processed
    // If pivot point (which should become a 1) is not yet on the diagonal, 
    // swap rows in the matrix equation (both matrix a and b). Note
    // that this is possible, since the order of the equation is
    // irrelevant. When swapping the rows, the solution doesn't change, so
    // the we do not need any book keeping to keep trackk of the changes we
    // made.
    if(pivot_row != pivot_col) {  
      ptr1 = matrixa + pivot_row*n; // points to the start of one of the rows to be swapped
      ptr2 = matrixa + pivot_col*n; // points to the start of the other row to be swapped
      for(i = 0; i < n; i++) {
	tmpvalue = *ptr1;
	*ptr1 = *ptr2;
	*ptr2 = tmpvalue;
      }
      // Now do the same for the right-hand side of the matrix equation
      ptr1 = matrixb + pivot_row*n; // points to the start of one of the rows to be swapped
      ptr2 = matrixb + pivot_col*n; // points to the start of the other row to be swapped
      for(i = 0; i < n; i++) {
	tmpvalue = *ptr1;
	*ptr1 = *ptr2;
	*ptr2 = tmpvalue;
      }
    }
    // So now the pivot element is arrived at the diagonal.


    // Rather than swapping columns, we keep track of what was implicitly swapped. 
    // The swapping will be done after the fact.
    orig_pivot_row[pivotnr] = pivot_row; 
    orig_pivot_column[pivotnr] = pivot_col;
    // Take the column in matrix a and scale the rows of a and b it so the pivot 
    // element becomes 1. Clearly dividing one of the set of equations we want
    // to solve (a row) with a constant does not change the solution.
    // To replace matrix a with the inverse matrix 
    // (rather than building it from the unity matrix),
    // we replace the diagonal element by 1 before we do the scaling.
    matrixa[pivot_col*n+pivot_col]=1.0;
    ptr1 = matrixa + pivot_col*n; // points to the start of the rows to be re-scaled
    for(i = 0; i < n; i++) {
      *ptr1 /= largest_element_value;
      ptr1 += 1; // Go to the next element
    }
    // Do the equivalent scaling with matrix b
    ptr1 = matrixb + pivot_col*m; // points to the start of the rows to be re-scaled
    for(i = 0; i < m; i++) {
      *ptr1 /= largest_element_value;
      ptr1 += 1; // Go to the next element
    }

    // We want to make the other elements of the column 0, while 
    // only the diagonal element is 1. This can be done by taking
    // the pivot row (which now has a 1 on the diagonal), and 
    // subtract it the the correct scaling from the other rows.
    // So loop over the rows.
    for(currow = 0; currow < n; currow++) {
      if(currow != pivot_col) {
	i = currow*n+pivot_col;
	tmpvalue = matrixa[i];  // This is the value that we want to make zero
	matrixa[i] = 0.0;       // Set this element as it would be in the identity matrix
	ptr1 = matrixa + n*currow;
	ptr2 = matrixa + n*pivot_col;
	for(i=1; i <= n; i++) {
	  *ptr1 -= (*ptr2)*tmpvalue;
	  ptr1++;
	  ptr2++;
	}
	ptr1 = matrixb + currow*m;
	ptr2 = matrixb + pivot_col*m;
	for(i = 0; i < m; i++) {
	  *ptr1 -= (*ptr2)*tmpvalue;
	  ptr1++;
	  ptr2++;
	}
      }
    }
  }  // End of the loop over the columns. So once done with the loop the solution is determined.

  // The column swaps were not applied yet, do the swapping now, 
  // since it is affect the order of the inverse of matrix a.
  for(coli = n-1; coli >= 0; coli--) {
    if(orig_pivot_row[coli] != orig_pivot_column[coli]) {
      for(rowi = 0; rowi <= n; rowi++) {
	tmpvalue = matrixa[rowi*n+orig_pivot_row[coli]];
	matrixa[rowi*n+orig_pivot_row[coli]] = matrixa[rowi*n+orig_pivot_column[coli]];
	matrixa[rowi*n+orig_pivot_column[coli]] = tmpvalue;
      }
    }
  }

  free(pivot_column_done);
  free(orig_pivot_row);
  free(orig_pivot_column);
  return 0;
}

//START REGION DEVELOP
