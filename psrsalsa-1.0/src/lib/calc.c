/*
calc_search_math_function assumes there are no letters in variable names

*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "psrsalsa.h"

#define calc_maxstringlength  1000
#define PlaceholderChar   '~'

double placeholders_internal[calc_maxstringlength];
int firstfreeplaceholder_internal;

int calc_internal_getval(char *equation, int start1, int end1, int nrvariables, char variables[][100], double *varables_values, double *value);
int calc_internal_getvals(char *equation, int start1, int end1, int start2, int end2, int nrvariables, char variables[][100], double *varables_values, double *value1, double *value2);
int calc_isanormalnumberchar(char c);
int calc_isanormalnumberorplaceholderchar(char c);
int calc_isanormalnumberstring(char *txt);
int calc_search_single_operator(char *equation, char operator, int *start1, int *end1, int *start2, int *end2);
int calc_search_matching_parenthesis(char *equation, int start, int *end);
int calc_search_parenthesis_open(char *equation, int *start);
void calc_insert_placeholder(char *equation, int start, int end);
int calc_search_math_function(char *equation, char *equationNew, char *function, int verbose);
int calc_substitute_variables_internal(char *valuestr, int nrvariables, char variables[][100], double *varables_values, double *value, int only_check_variable);

/* Print the mathematics functions to device (could be for instance
   stdio). nrspaces defines the number of spaces before each line.*/
void printCalcFunctions(FILE *printdevice, int nrspaces)
{
  int i;
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "Opperators: +,-,*,/,^\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "sin\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "cos\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "tan\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "asin\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "acos\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "atan\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "sinh\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "cosh\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "tanh\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "sqrt\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "exp\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "log or ln\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "log10\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "deg          (make a value between 0 and 360)\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "degsmall     (179->179, 181->-179, 270->-90)\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "deg180       (make a value between 0 and 180)\n");
}

int calc_expression_internal(char *equation, int loopnr, int nrvariables, char variables[][100], double *varables_values, int verbose)
{
  int ret, done, i, startpar, endpar, start1, end1, start2, end2, tmpplaceholder;
  double value1, value2, value3;
  char newequation[calc_maxstringlength], function[calc_maxstringlength];

  done = 0;
  if(verbose) {
    if(loopnr == 0) {
      printwarning(0, "\ncalc_expression_internal (start of new call itt=%d):  %s", loopnr+1, equation); 
    }else {
      printf("calc_expression_internal (start of new call itt=%d):  %s\n", loopnr+1, equation); 
    }
  }
  do {
    if(verbose) printf("calc_expression_internal (next loop of itt=%d):  %s\n", loopnr+1, equation); 
    /* Search for functions first  */
    if(calc_search_math_function(equation, newequation, function, verbose)) {
      if(verbose) printf("  found function '%s' (%s)\n", function, newequation); 
      /* Remember the placeholder that will be used to store the output value of the function */
      tmpplaceholder = firstfreeplaceholder_internal-1;
      if(verbose) printf("  Request calculation of function argument '%s' (itt=%d)\n", newequation, loopnr+1);
      /* Calculate the input value of the function by calling calc_expression_internal recursively */
      ret = calc_expression_internal(newequation, loopnr+1, nrvariables, variables, varables_values, verbose);
      if(ret == 0) {
	fflush(stdout);
	fprintf(stderr, "ERROR calc_expression_internal: Failed to evaluate '%s'.\n", newequation);
	return 0;
      }
      /* The last computation done must be the answer of the input expression */
      value1 = placeholders_internal[firstfreeplaceholder_internal-1];
      if(verbose) printf("  Answer of function argument calculation '%s': %e (%c%d)    (itt=%d)\n", newequation, value1, PlaceholderChar, firstfreeplaceholder_internal-1, loopnr+1);
      /* Compute the output value */
      if(strcmp(function, "sin") == 0) {
	value2 = sin(value1);
      }else if(strcmp(function, "cos") == 0) {
	value2 = cos(value1);
      }else if(strcmp(function, "tan") == 0) {
	value2 = tan(value1);
      }else if(strcmp(function, "asin") == 0) {
	value2 = asin(value1);
      }else if(strcmp(function, "acos") == 0) {
	value2 = acos(value1);
      }else if(strcmp(function, "atan") == 0) {
	value2 = atan(value1);
      }else if(strcmp(function, "sinh") == 0) {
	value2 = sinh(value1);
      }else if(strcmp(function, "cosh") == 0) {
	value2 = cosh(value1);
      }else if(strcmp(function, "tanh") == 0) {
	value2 = tanh(value1);
	/*
      }else if(strcmp(function, "asinh") == 0) {
	value2 = asinh(value1);
      }else if(strcmp(function, "acosh") == 0) {
	value2 = acosh(value1);
      }else if(strcmp(function, "atanh") == 0) {
	value2 = atanh(value1);
	*/
      }else if(strcmp(function, "sqrt") == 0) {
	value2 = sqrt(value1);
      }else if(strcmp(function, "exp") == 0) {
	value2 = exp(value1);
      }else if(strcmp(function, "log") == 0 || strcmp(function, "ln") == 0) {
	value2 = log(value1);
      }else if(strcmp(function, "log10") == 0) {
	value2 = log10(value1);
      }else if(strcmp(function, "deg") == 0) {
	value2 = derotate_deg_double(value1); 
      }else if(strcmp(function, "degsmall") == 0) {
	value2 = derotate_deg_small_double(value1);
      }else if(strcmp(function, "deg180") == 0) {
	value2 = derotate_180_double(value1);
     }else {
	fflush(stdout);
	fprintf(stderr, "ERROR calc_expression_internal: Don't recognize function '%s'.\n", function);
	return 0;
      }
      placeholders_internal[tmpplaceholder] = value2;
      if(verbose) printf("  %s(%e)=%e (stored as %c%d)\n", function, value1, value2, PlaceholderChar, tmpplaceholder);
      //      if(verbose) printf("  %s(%f)=%f\n", function, value1, value2);
    }else if(calc_search_parenthesis_open(equation, &startpar)) {  /* Are there expressions in parenthesis? */
      if(!calc_search_matching_parenthesis(equation, startpar, &endpar)) {
	fflush(stdout);
	fprintf(stderr, "ERROR calc_expression_internal: Cannot find closing ).\n");
	return 0;
      }
      strcpy(newequation, &equation[startpar+1]);
      newequation[endpar-startpar-1] = 0;
      /* This is the placeholder that will be used for the result in parenthesis */
      tmpplaceholder = firstfreeplaceholder_internal;
      /* Replace the expression in parenthesis by a placehoder */
      calc_insert_placeholder(equation, startpar, endpar);
      if(verbose) printf("  going to calculate '%s' (%s) as in parenthesis  (itt=%d)\n", newequation, equation, loopnr+1);
      //      if(verbose) printf("  going to calculate '%s'\n", newequation);
      /* Recursively calculate expression in parenthesis */
      ret = calc_expression_internal(newequation, loopnr+1, nrvariables, variables, varables_values, verbose);
      if(ret == 0) {
	fflush(stdout);
	fprintf(stderr, "ERROR calc_expression_internal: Failed to evaluate '%s'.\n", newequation);
	return 0;
      }
      /* The last computation must be the answer of the expression in parenthesis */
      value1 = placeholders_internal[firstfreeplaceholder_internal-1];
      placeholders_internal[tmpplaceholder] = value1;
      if(verbose) printf("Put %e in %c%d from %c%d (expression in parenthesis, itt=%d)\n", value1, PlaceholderChar, tmpplaceholder, PlaceholderChar, firstfreeplaceholder_internal-1, loopnr+1);
    }else if(calc_search_single_operator(equation, '^', &start1, &end1, &start2, &end2)) {  
      if(verbose) {
	printf("  Found ^ operator acting on ");
	printf("'");
	for(i = start1; i <= end1; i++) {
	  printf("%c", equation[i]);
	}
	printf("' and '");
	for(i = start2; i <= end2; i++) {
	  printf("%c", equation[i]);
	}
	printf("'\n");
      }
      if(calc_internal_getvals(equation, start1, end1, start2, end2, nrvariables, variables, varables_values, &value1, &value2) == 0) {
	fflush(stdout);
	fprintf(stderr, "ERROR calc_expression_internal: Failed to evaluate '%s'.\n", equation);
	return 0;
      }
      calc_insert_placeholder(equation, start1, end2);
      value3 = pow(value1, value2);
      if(verbose) printf("  Calculating '%f^%f=%f' (%c%d)\n", value1, value2, value3, PlaceholderChar, firstfreeplaceholder_internal-1);
      //      if(verbose) printf("  Calculating '%f^%f = %f'\n", value1, value2, value3); 
      placeholders_internal[firstfreeplaceholder_internal-1] = value3;
    }else if(calc_search_single_operator(equation, '*', &start1, &end1, &start2, &end2)) {  
      if(verbose) {
	printf("  Found * operator acting on ");
	printf("'");
	for(i = start1; i <= end1; i++) {
	  printf("%c", equation[i]);
	}
	printf("' and '");
	for(i = start2; i <= end2; i++) {
	  printf("%c", equation[i]);
	}
	printf("'\n");
      }
      if(calc_internal_getvals(equation, start1, end1, start2, end2, nrvariables, variables, varables_values, &value1, &value2) == 0) {
	fflush(stdout);
	fprintf(stderr, "ERROR calc_expression_internal: Failed to evaluate '%s'.\n", equation);
	return 0;
      }
      calc_insert_placeholder(equation, start1, end2);
      value3 = value1 * value2;
      if(verbose) printf("  Calculating '%f*%f=%f' (%c%d)\n", value1, value2, value3, PlaceholderChar, firstfreeplaceholder_internal-1); 
      //      if(verbose) printf("  Calculating '%f*%f = %f'\n", value1, value2, value3); 
      placeholders_internal[firstfreeplaceholder_internal-1] = value3;
    }else if(calc_search_single_operator(equation, '/', &start1, &end1, &start2, &end2)) {  
      if(verbose) {
	printf("  Found / operator acting on ");
	printf("'");
	for(i = start1; i <= end1; i++) {
	  printf("%c", equation[i]);
	}
	printf("' and '");
	for(i = start2; i <= end2; i++) {
	  printf("%c", equation[i]);
	}
	printf("'\n");
      }
      if(calc_internal_getvals(equation, start1, end1, start2, end2, nrvariables, variables, varables_values, &value1, &value2) == 0) {
	fflush(stdout);
	fprintf(stderr, "ERROR calc_expression_internal: Failed to evaluate '%s'.\n", equation);
	return 0;
      }
      calc_insert_placeholder(equation, start1, end2);
      value3 = value1 / value2;
      if(verbose) printf("  Calculating '%f/%f=%f' (%c%d)\n", value1, value2, value3, PlaceholderChar, firstfreeplaceholder_internal-1);
      //      if(verbose) printf("  Calculating '%f/%f = %f'\n", value1, value2, value3); 
      placeholders_internal[firstfreeplaceholder_internal-1] = value3;
      /* Do minus first, or else 1-2+3 = -4 */
    }else if(calc_search_single_operator(equation, '-', &start1, &end1, &start2, &end2)) {  
      if(end1 >= 0) {
	if(verbose) {
	  printf("  Found - operator acting on ");
	  printf("'");
	  for(i = start1; i <= end1; i++) {
	    printf("%c", equation[i]);
	  }
	  printf("' and '");
	  for(i = start2; i <= end2; i++) {
	    printf("%c", equation[i]);
	  }
	  printf("'\n");
	}
	if(calc_internal_getvals(equation, start1, end1, start2, end2, nrvariables, variables, varables_values, &value1, &value2) == 0) {
	  fflush(stdout);
	  fprintf(stderr, "ERROR calc_expression_internal: Failed to evaluate '%s'.\n", equation);
	  return 0;
	}
	calc_insert_placeholder(equation, start1, end2);
	value3 = value1 - value2;
	if(verbose) printf("  Calculating '%f-%f=%f' (%c%d)\n", value1, value2, value3, PlaceholderChar, firstfreeplaceholder_internal-1);
	//	if(verbose) printf("  Calculating '%f-%f = %f'\n", value1, value2, value3); 
	placeholders_internal[firstfreeplaceholder_internal-1] = value3;
      }else {
	/* Only one value like -1 (instead of something like 3-1). So get value and multimply by -1. */
	if(verbose) {
	  printf("  Found - (making neg number) operator acting on ");
	  printf("'");
	  for(i = start2; i <= end2; i++) {
	    printf("%c", equation[i]);
	  }
	  printf("'\n");
	}
	if(calc_internal_getval(equation, start2, end2, nrvariables, variables, varables_values, &value1) == 0) {
	  fflush(stdout);
	  fprintf(stderr, "ERROR calc_expression_internal: Failed to evaluate '%s'.\n", equation);
	  return 0;
	}
	calc_insert_placeholder(equation, start1, end2);
	value1 *= -1.0;
	if(verbose) printf("  Calculating '%f = %f' (%c%d)\n", value1, value1, PlaceholderChar, firstfreeplaceholder_internal-1);
	//	if(verbose) printf("  Calculating '%f = %f'\n", value1, value1); 
	placeholders_internal[firstfreeplaceholder_internal-1] = value1;
      }
    }else if(calc_search_single_operator(equation, '+', &start1, &end1, &start2, &end2)) {  
      if(verbose) {
	printf("  Found + operator acting on ");
	printf("'");
	for(i = start1; i <= end1; i++) {
	  printf("%c", equation[i]);
	}
	printf("' and '");
	for(i = start2; i <= end2; i++) {
	  printf("%c", equation[i]);
	}
	printf("'\n");
      }
      if(calc_internal_getvals(equation, start1, end1, start2, end2, nrvariables, variables, varables_values, &value1, &value2) == 0) {
	fflush(stdout);
	fprintf(stderr, "ERROR calc_expression_internal: Failed to evaluate '%s'.\n", equation);
	return 0;
      }
      calc_insert_placeholder(equation, start1, end2);
      value3 = value1 + value2;
      if(verbose) printf("  Calculating '%f+%f=%f' (%c%d)\n", value1, value2, value3, PlaceholderChar, firstfreeplaceholder_internal-1);
      //      if(verbose) printf("  Calculating '%f+%f = %f'\n", value1, value2, value3); 
      placeholders_internal[firstfreeplaceholder_internal-1] = value3;
    }else if(calc_isanormalnumberstring(equation)) {  /* Is it just a number? */
      sscanf(equation, "%lf", &value1);
      calc_insert_placeholder(equation, 0, strlen(equation));
      placeholders_internal[firstfreeplaceholder_internal-1] = value1;
      if(verbose) printf("Put number '%e' (%c%d)\n", value1, PlaceholderChar, firstfreeplaceholder_internal-1);
    }else if(calc_substitute_variables_internal(equation, nrvariables, variables, varables_values, &value1, 1)){ /*  is it a variable? */
      done = 1;
      /* This is the placeholder that will be used for the result in parenthesis */
      tmpplaceholder = firstfreeplaceholder_internal;
      /* Replace the expression in parenthesis by a placehoder */
      calc_insert_placeholder(equation, 0, strlen(equation)-1);
      placeholders_internal[tmpplaceholder] = value1;
      if(verbose) printf("Put '%e' (%c%d)\n", value1, PlaceholderChar, tmpplaceholder);
    }else if(equation[0] == PlaceholderChar) {  /* Is it just a placeholder? Then the expression is calculated. */
      done = 1;
      /* Do a quick test if it is really just a placeholder, if not, there is a bug. */
      for(i = 1; i < strlen(equation); i++) {
	if(equation[i] < '0' || equation[i] > '9') {
	  fflush(stdout);
	  fprintf(stderr, "ERROR calc_expression_internal: Bug! Expected a placeholder '%s'.\n", equation);
	  return 0;
	}
      }
      if(verbose) printf("Just a placeholder, must be done (itt=%d; eq='%s')\n", loopnr+1, equation);
      int value_placeholder;
      sscanf(&equation[1], "%d", &value_placeholder);
      // If not doing this, the calculation of an expression like: (sin($1*3.1415/180)) goes wrong, as answer of the argument of the outer () is not stored in the last spot, while it is assumed that the result of the argument of () is in the last spot.
      if(value_placeholder != firstfreeplaceholder_internal - 1) {
	if(verbose) printf("  but placeholder is not last placeholder. Make a copy\n");
	if(verbose) printf("  Put %e in %c%d\n", placeholders_internal[value_placeholder], PlaceholderChar, firstfreeplaceholder_internal);
	placeholders_internal[firstfreeplaceholder_internal] = placeholders_internal[value_placeholder];
	firstfreeplaceholder_internal++;
      }
    }else { /* Something is wrong with expression, unless it is a variable. */
      fflush(stdout);
      fprintf(stderr, "ERROR calc_expression_internal: syntax error '%s'.\n", equation);
      return 0;
    }
  }while(done == 0);   /* Keep computing until solved. */
  return 1;
}

/* Calculates equation and puts it in answer. Set verbose to 1 to see
   the answer and set it to 2 to show all intermediate steps. All
   calculations are done in double format, but the answer can both be
   a float or a double. Returns 0 on error. */
int calc_expression(char *equation, int nrvariables, char variables[][100], double *varables_values, double *answer, int verbose)
{
  char newequation[calc_maxstringlength];
  int ret, verbose2, i, i2;
  double value;
  i2 = 0;
  for(i = 0; i < strlen(equation); i++) {
    if(equation[i] != ' ' && equation[i] != '\t' && equation[i] != '\n' && equation[i] != '\r') {
      newequation[i2++] = equation[i];
    }
  }
  newequation[i2] = 0;
  firstfreeplaceholder_internal = 0;
  if(verbose == 1 || verbose == 0)
    verbose2 = 0;
  else
    verbose2 = 1;
  ret = calc_expression_internal(newequation, 0, nrvariables, variables, varables_values, verbose2);
  if(verbose2)
    printf("final expression stored as='%s'\n", newequation);
  ret += calc_internal_getval(newequation, 0, strlen(newequation)-1, nrvariables, variables, varables_values, &value);
  *answer = value;
  if(ret == 2) {
    if(verbose) printf("Final expression: '%s=%f'\n", equation, value);
    return 1;
  }else {
    if(verbose) printf("ERROR calc_expression: Failed to calculate '%s'\n", equation);
    return 0;
  }
  return ret;
}

int calc_expressionf(char *equation, int nrvariables, char variables[][100], float *varables_values, float *answerf, int verbose)
{
  int ret;
  double *values, answer;
  if(nrvariables > 0) {
    values = (double *)malloc(nrvariables*sizeof(double));
    if(values == NULL) {
      fflush(stdout);
      fprintf(stderr, "ERROR calc_expressionf: Memory allocation error.\n");
      return 0;
    }
  }
  for(ret = 0; ret < nrvariables; ret++)
    values[ret] = varables_values[ret];
  ret = calc_expression(equation, nrvariables, variables, values, &answer, verbose);
  *answerf = answer;
  if(nrvariables > 0) {
    free(values);
  }
  return ret;
}
/*
int main()
{
  char txt[calc_maxstringlength], *txtptr;
  double answer;
  do {
    txtptr = fgets(txt, calc_maxstringlength, stdin);
    if(txtptr != NULL) {
      if(txtptr[strlen(txtptr)-1] == '\n')
	txtptr[strlen(txtptr)-1] = 0;
      calc_expression(txt, &answer, 2);
    }
  }while(txtptr != NULL);
  return 0;
}
*/

int calc_substitute_variables_internal(char *valuestr, int nrvariables, char variables[][100], double *varables_values, double *value, int only_check_variable)
{
  int i, ret;
  ret = 0;
  // See if value really is a known variable
  for(i = 0; i < nrvariables; i++) {
    if(strcmp(valuestr, variables[i]) == 0) {
      *value = varables_values[i];
      /*      printf("subs: '%s'=%f\n", valuestr, *value); */
      ret = 1;
    }
  }
  // See if variable is a build in variable, like pi=3.1415...
  if(strcmp(valuestr, "pi") == 0 || strcmp(valuestr, "Pi") == 0 || strcmp(valuestr, "PI") == 0) {
    *value = M_PI;
    ret = 1;
  }
  if(only_check_variable)
    return ret;

  // No, then it should be a value in ascii
  if(ret == 0) {
    if(calc_isanormalnumberstring(valuestr) == 0) {
      printerror(0, "ERROR calc_substitute_variables_internal: Cannot interpret '%s' as a value.", valuestr);
      return 0;
    }
    ret = sscanf(valuestr, "%lf", value);
  }
  if(ret != 1) {
    /*      printf("Oeps: '%s'\n", valuestr); */
    return 0;
  }
  return 1;
}

int calc_internal_getval(char *equation, int start1, int end1, int nrvariables, char variables[][100], double *varables_values, double *value)
{
  int ret, index;
  char valuestr[calc_maxstringlength];
  strcpy(valuestr, &equation[start1]);
  valuestr[end1-start1+1] = 0;
  /*  printf("%s %d %d\n", valuestr, start1, end1); */
  if(valuestr[0] == PlaceholderChar) {
    // The placeholder char '~' should be followed by just numbers only (integer), nothing else (spaces are already removed)
    if(strlen(valuestr) < 2) {
      printerror(0, "ERROR calc_internal_getval: '%s' cannot be interpreted as a placeholder.", valuestr);
      return 0;
    }
    for(index = 1; index < strlen(valuestr); index++) {
      if(valuestr[index] < '0' || valuestr[index] > '9') {
	printerror(0, "ERROR calc_internal_getval: '%s' cannot be interpreted as a placeholder.", valuestr);
	return 0;
      }
    }
    ret = sscanf(&valuestr[1], "%d", &index);
    if(ret != 1)
      return 0;
    *value = placeholders_internal[index];
  }else {
    ret = calc_substitute_variables_internal(valuestr, nrvariables, variables, varables_values, value, 0);
    if(ret == 0)
      return 0;
  }
  return 1;
}




int calc_internal_getvals(char *equation, int start1, int end1, int start2, int end2, int nrvariables, char variables[][100], double *varables_values, double *value1, double *value2)
{
  int ret;
  ret  = calc_internal_getval(equation, start1, end1, nrvariables, variables, varables_values, value1);
  ret += calc_internal_getval(equation, start2, end2, nrvariables, variables, varables_values, value2);
  if(ret != 2)
    return 0;
  return 1;
}



int calc_isanormalnumberchar(char c)
{
  if((c >= '0' && c <= '9') || c == '.')
    return 1;
  else
    return 0;
}

int calc_isanormalnumberorplaceholderchar(char c)
{
  if(calc_isanormalnumberchar(c)) {
    return 1;
  }else {
    if(c == PlaceholderChar)
      return 1;
    else
      return 0;
  }
}

int calc_isanormalnumberstring(char *txt)
{
  int i;
  int number_e, number_dot, number_plusmin;
  number_e = 0;
  number_dot = 0;
  number_plusmin = 0;
  for(i = 0; i < strlen(txt); i++) {
    if(txt[i] < '0' || txt[i] > '9') { // Not a number
      if(txt[i] == 'e' || txt[i] == 'E') {
	number_e += 1;
      }else if(txt[i] == '.') {
	number_dot += 1;
      }else if(txt[i] == '+' || txt[i] == '-') {
	number_plusmin += 1;
      }else {
	return 0;
      }
    }
  }
  if(number_dot > 1) {
    return 0;
  }
  if(number_e > 1) {
    return 0;
  }
  if(number_e == 1) {
    if(number_plusmin > 1) {
	return 0;
    }
  }else {
    if(number_plusmin > 0) {
      return 0;
    }
  }
  //  for(i = 0; i < strlen(txt); i++) 
  //    if(!calc_isanormalnumberchar(txt[i]))
  //      return 0;
  return 1;
}


// Check if + or - (at location equation[plusminus_pos]) is part of exponential (1.23e-4)
// return 1=yes, 0 = no
int check_plusminus_part_exponential(char *equation, int plusminus_pos)
{
  // There must be a character following
  if(plusminus_pos == strlen(equation) - 1) {
    return 0;
  }
  // Must be preceeded by at least an e and a number
  if(plusminus_pos < 2) {
    return 0;
  }
  // Must be followed by a number
  if(equation[plusminus_pos+1] < '0' || equation[plusminus_pos+1] > '9') {
    return 0;
  }
  // Must be preceeded by 'e'
  if(equation[plusminus_pos-1] != 'e' && equation[plusminus_pos-1] != 'E') {
    return 0;
  }
  // Must be preceeded by a number or a .
  if((equation[plusminus_pos-2] < '0' || equation[plusminus_pos-2] > '9') && equation[plusminus_pos-2] != '.') {
    return 0;
  }
  return 1;
}

int calc_search_single_operator(char *equation, char operator, int *start1, int *end1, int *start2, int *end2)
{
  int found = 0, ioperator;
  // Look if operator is present in string.
  for(ioperator = 0; ioperator < strlen(equation); ioperator++) {
    if(equation[ioperator] == operator) {
      if(operator == '+' || operator == '-') {  // Make sure +/- not part of exponential, i.e. 1.23e+4
	if(check_plusminus_part_exponential(equation, ioperator) == 0) {
	  found = 1;
	  break;
	}
      }else {
	found = 1;
	break;
      }
    }
  }
  // Now isolate a substring from the equation of the form ..... operator ......
  if(found) {
    /*    printf("Found '%c' in '%s' at %d\n", operator, equation, ioperator); */
    // Find the start of the substring
    // Note: (,) already taken care of before looking for operators, so there shouldn't be any in the string.
    // Same is true for functions, so it should be actual numbers + placeholders + variables.
    // Since variable names are unknown, look for other operators which act as boundaries
    *end1 = ioperator - 1;
    for(*start1 = ioperator-1; *start1 >= 0; (*start1)--) {
      /*
      if(!((equation[*start1] >= '0' && equation[*start1] <= '9') || equation[*start1] == PlaceholderChar || equation[*start1] == '.')) {
	break;
      }
      */
      if(equation[*start1] == '^' || equation[*start1] == '*' || equation[*start1] == '/') {
	break;
      }
      //These could be boundaries, or part of exponentials
      if(equation[*start1] == '+' || equation[*start1] == '-') {
	if(check_plusminus_part_exponential(equation, *start1) == 0) // Not part of exponential, then it is a boundary
	  break;
      }
    }
    (*start1) += 1;
    *start2 = ioperator+1;
    for(*end2 = ioperator+1; *end2 <= strlen(equation); (*end2)++) {
      /*      if(!((equation[*end2] >= '0' && equation[*end2] <= '9') || equation[*end2] == PlaceholderChar || equation[*end2] == '.')) {
	break;
	}*/
      if(equation[*end2] == '^' || equation[*end2] == '*' || equation[*end2] == '/'  || equation[*end2] == 0) {
	break;
      }
      if(equation[*end2] == '+' || equation[*end2] == '-') {
	if(check_plusminus_part_exponential(equation, *end2) == 0) // Not part of exponential, then it is a boundary
	  break;
      }
    }
    (*end2) -= 1;
    /*    printf("%d %d\n", *start1, *end2); */
  }
  return found;
}

int calc_search_matching_parenthesis(char *equation, int start, int *end)
{
  int found = 0, depth;
  if(equation[start] != '(') {
    fflush(stdout);
    fprintf(stderr, "ERROR calc_search_matching_parenthesis: Doesn't have a ( at expected location.\n");
  }else {
    depth = 1;
    for(*end = start+1; *end < strlen(equation); (*end)++) {
      if(equation[*end] == '(') {
	depth++;
      }else if(equation[*end] == ')') {
	depth--;
	if(depth == 0) {
	  found = 1;
	  break;
	}
      }
    }
  }
  return found;
}

int calc_search_parenthesis_open(char *equation, int *start)
{
  int found = 0;
  for(*start = 0; *start < strlen(equation); (*start)++) {
    if(equation[*start] == '(') {
      found = 1;
      break;
    }
  }
  return found;
}

void calc_insert_placeholder(char *equation, int start, int end)
{
  char txt[100], tmpequation[calc_maxstringlength];
  int i, j, k;

  /*  printf("HOI: %d %d '%s'\n", start, end, equation); */
  strcpy(tmpequation, equation);
  k = strlen(equation);
  equation[start] = PlaceholderChar;
  equation[start+1] = 0;
  sprintf(txt, "%d", firstfreeplaceholder_internal++);
  strcat(equation, txt);
  /*  printf("HOIHOI %d %d\n", end+1, (int)strlen(equation)); */
  for(i = end+1; i < k; i++) {
    j = strlen(equation);
    equation[j] = tmpequation[i];
    equation[j+1] = 0;
  }
}

int calc_search_math_function(char *equation, char *equationNew, char *function, int verbose)
{
  int startfunc, endfunc, starteqnew, endeqnew, found;
  found = 0;
  for(startfunc = 0; startfunc < strlen(equation); startfunc++) {
    // Look for alphabet character which could be the start of a function
    if((equation[startfunc] >= 'a' && equation[startfunc] <= 'z') || (equation[startfunc] >= 'A' && equation[startfunc] <= 'Z')) {
      for(endfunc = startfunc+1; endfunc < strlen(equation); endfunc++) {
	// See if first alphabet character is followed by maybe more characters and there should be a parenthesis open. Numerical characters are also allowed (i.e. log10).
	if(equation[endfunc] == '(') {
	  found = 1;
	  break;
	}else if((equation[endfunc] >= 'a' && equation[endfunc] <= 'z') || (equation[endfunc] >= 'A' && equation[endfunc] <= 'Z') || (equation[endfunc] >= '0' && equation[endfunc] <= '9')) {
	  /* These characters are ok in function name */
	  //	  printf("XXXX accepted '%c' in function name\n", equation[endfunc]);
	}else {
	  /* Not a function name */
	  found = 0;
	  break;
	}
      }
    }
    if(found)
      break;
  }
  if(found) {
    strncpy(function, &equation[startfunc], endfunc-startfunc);
    function[endfunc-startfunc] = 0;
    starteqnew = endfunc+1;
    if(calc_search_matching_parenthesis(equation, endfunc, &endeqnew)) {
      strcpy(equationNew, &equation[starteqnew]);
      equationNew[endeqnew-starteqnew] = 0;
      calc_insert_placeholder(equation, startfunc, endeqnew);
      /*
	printf("%d %d %d %d\n", startfunc, endfunc, starteqnew, endeqnew);
	printf("%s\n", equationNew);
	printf("%s\n", equation);
      */
    }else {
      fflush(stdout);
      fprintf(stderr, "ERROR calc_search_math_function: Function '%s' doesn't have a )", function);
    }
  }
  return found;
}
