#include <stdio.h>
#include <R.h>
#include <Rdefines.h>
#include <vector>
#include <list>
#include <bitset>
#include <cmath>

typedef std::list< std::vector<double> > vector_list_t;

// From: http://stackoverflow.com/questions/16761472/how-can-i-increment-stdbitset
template <size_t N>
void increment(std::bitset<N>& in)
{
	for ( size_t i = 0; i < N; ++i )
	{
		if ( in[i] == 0 )
		{
			in[i] = 1;
			break;
		}
		in[i] = 0;
	}
}

void find_vertices(vector_list_t& vert, double p, const double* Pi, int J)
{
	int r = 0;
	int bin_cnt = (int) pow(2.0, (double)(J-1) );
	int q = 0;

	for (int j = 0; j < J; j++)
	{
		// If Pi[j] == 0, we can't divide by it and get a solution for v.j 
		if (Pi[j] == 0)
			continue;

		// x <- as.vector( digitsBase(k, base = 2, J-1) )
		for (std::bitset<64> x(0); x.to_ulong() < bin_cnt; increment(x))
		{
			double ss = 0;
			q = 0;
			for (int l = 0; l < J; l++)
			{
				if (l != j)
				{
					// If the q-th digit (base 2) of x is == 1, add Pi[l] to ss
					ss += Pi[l] * x[q];
					q++;
				}
			}

			// v.j <- 1/Pi[j] * (p - t(x) %*% Pi[-j])
			double v_j = 1/Pi[j] * (p - ss);

			if (0 <= v_j && v_j <= 1)
			{
				r++;
				std::vector<double> v(J);
			
				// v[-j] <- x
				// v[j] <- v.j

				q = 0;
				for (int l = 0; l < J; l++)
				{
					if (l == j)
					{
						v[l] = v_j;
					}
					else
					{
						v[l] = x[q];
						q++;
					}
				}

				vert.push_back(v);
			}
		}
	}
}

double norm2(const std::vector<double>& x, const std::vector<double>& y)
{
	double ss = 0;
	int n = x.size();
	int n2 = y.size();

	// TBD: Throw an error if n != n2

	for (int i = 0; i < n; i++)
	{
		ss += pow(x[i] - y[i], 2.0);
	}

	return sqrt(ss);
}

void extract_unique_vertices(vector_list_t& unique_vert_list,
	const vector_list_t& vert_list, double tol)
{
	for (vector_list_t::const_iterator itr = vert_list.begin();
		itr != vert_list.end(); itr++)
	{
		// Compare the vertex at itr to all previous vertices
		// in the unique list. If norm distance is very small, count it
		// as a match
		bool is_dup = FALSE;

		for (vector_list_t::const_iterator uniq_itr = unique_vert_list.begin();
			uniq_itr != unique_vert_list.end(); uniq_itr++)
		{
			double d = norm2(*itr, *uniq_itr);
			if (d < tol)
			{
				is_dup = TRUE;
				break;
			}
		}

		if (!is_dup)
		{
			unique_vert_list.push_back(*itr);
		}
	}
}

void vector_list_to_matrix(const vector_list_t& vert_list,
	double* V, int J, int k)
{
	int l = 0;
	for (vector_list_t::const_iterator itr = vert_list.begin();
		itr != vert_list.end(); itr++)
	{
		const std::vector<double>& v = *itr;
		for (int j = 0; j < J; j++)
		{
			V[l*J + j] = v[j];
		}

		l++;
	}
}

// This function just packs and unpacks data for C++
extern "C"
SEXP find_vertices(SEXP sexp_p, SEXP sexp_Pi, SEXP sexp_tol)
{
	double p = NUMERIC_VALUE(sexp_p);
	double* Pi = NUMERIC_POINTER(sexp_Pi);
	int J = GET_LENGTH(sexp_Pi);
	double tol = NUMERIC_VALUE(sexp_tol);

	// Call the C++ function to do the real work
	vector_list_t vert_list;
	vector_list_t unique_vert_list;

	find_vertices(vert_list, p, Pi, J);
	extract_unique_vertices(unique_vert_list, vert_list, tol);
	int k = unique_vert_list.size();

	// Now pack up the list of vertices into a J x k matrix and return it to R
	SEXP sexp_vert;
	PROTECT(sexp_vert = NEW_NUMERIC(J * k));
	double* V = NUMERIC_POINTER(sexp_vert);
	vector_list_to_matrix(unique_vert_list, V, J, k);

	SEXP sexp_dim;
	PROTECT(sexp_dim = NEW_INTEGER(2));
	int* dim = INTEGER_POINTER(sexp_dim);
	dim[0] = J;
	dim[1] = k;

	SET_DIM(sexp_vert, sexp_dim);

	UNPROTECT(2);
	return sexp_vert;
}

