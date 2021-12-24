#include "./headers/propagate_intervals.h"
#include "../flowstar-release/Continuous.h"

using namespace std;
using namespace flowstar;

datatype offset_in_constraint_comb = constr_comb_offset;

int main()
{
	Variables stateVars;

	/*
	 * Declaration of the state variables.
	 * The first one should always be the local time variable which can be
	 * viewed as a preserved variable. It is only used internally by the library.
	 */
	stateVars.declareVar("t");
	stateVars.declareVar("x");
	stateVars.declareVar("y");
	stateVars.declareVar("u");

	int domainDim = 4;

	Expression_AST deriv_x("y", stateVars);
	Expression_AST deriv_y("2 * sin(x) + 8 * u", stateVars);
	Expression_AST deriv_u("0", stateVars);

	ODE plant(stateVars);
	plant.assignDerivative("x", deriv_x);
	plant.assignDerivative("y", deriv_y);
	plant.assignDerivative("u", deriv_u);


	/*
	 * Specify the parameters for reachability computation.
	 */
	Continuous_Reachability_Setting crs;

	// step size
	crs.setFixedStepsize(0.01);

	// Taylor model order
	crs.setFixedOrder(4);

	// precision
	crs.setPrecision(100);

	// cutoff threshold
	Interval cutoff(-1e-10,1e-10);
	crs.setCutoff(cutoff);

	/*
	 * A remainder estimation is a vector of intervals such that
	 * the i-th component is the estimation for the i-th state variable.
	 */
	Interval E(-0.01,0.01);
	std::vector<Interval> estimation;
	estimation.push_back(E);	// estimation for the 1st variable
	estimation.push_back(E);	// estimation for the 2nd variable
	estimation.push_back(E);	// estimation for the 3rd variable
	crs.setRemainderEstimation(estimation);

	// call this function whenever a parameter is set or changed
	crs.prepareForReachability();

	// Simple range propagation
	char controller_file[] = "controller_Single-Pendulum_sherlock";

	network_handler system_network(controller_file);

	/*
	 * Initial set can be a box which is represented by a vector of intervals.
	 * The i-th component denotes the initial set of the i-th state variable.
	 */
	Interval init_x(1.0,1.2), init_y(0.0,0.2), init_u, intZero;
	std::vector<Interval> X0;
	X0.push_back(init_x);
	X0.push_back(init_y);
	X0.push_back(init_u);


	// translate the initial set to a flowpipe
	Flowpipe initial_set(X0, intZero);

	// the flowpipe that keeps the overapproximation at the end of a time horizon
	Flowpipe fp_last;

	// the symbolic remainder
	Symbolic_Remainder symb_rem(initial_set);

	std::list<Flowpipe> result;

	datatype max_difference = 0;

	std :: chrono :: steady_clock::time_point t1 = std :: chrono :: steady_clock::now();


	std:: chrono::duration< double > time_span;
	std :: chrono :: steady_clock::time_point start_time;
	std :: chrono :: steady_clock::time_point end_time;

	// std:: chrono::duration< double > total_time_spent_in_regression  = 0;
	double total_time_spent_in_regression  = 0;
	double total_time_spent_in_PWL_construction  = 0;
	double total_time_spent_in_Sherlock  = 0;
	double total_time_spent_in_Flowstar  = 0;
	int no_of_piecewise_linear_pieces= 0 ;
	int linear_piece_count = 0;
	int max_linear_pieces = 0;

	for(int k=0; k<20; ++k)
	{
		std::vector<Interval> NN_input;
		initial_set.intEvalNormal(NN_input, crs.step_end_exp_table, crs.cutoff_threshold);

		vector< vector< datatype > > input_interval(2, vector< datatype >(2,0));
		input_interval[0][0] = NN_input[0].inf();
		input_interval[0][1] = NN_input[0].sup();

		input_interval[1][0] = NN_input[1].inf();
		input_interval[1][1] = NN_input[1].sup();

		printf("[%lf, %lf], \t [%lf, %lf]\n",
		       input_interval[0][0], input_interval[0][1],
		       input_interval[1][0], input_interval[1][1]);

		vector< vector< datatype > > input_constraints;
		create_constraint_from_interval(input_constraints, input_interval);

		vector< vector< unsigned int > > monomial_terms;
		vector< datatype > coefficients;
		datatype offset = 0.0;
		datatype scaling = 1.0;
		unsigned int degree = 2;

		start_time = std :: chrono :: steady_clock::now();
		generate_polynomial_for_NN(system_network, degree, input_constraints,
					   offset, scaling, monomial_terms, coefficients);
		int i,j,l;
		vector<my_monomial_t> polynomial;
		polynomial = create_polynomial_from_monomials_and_coeffs(monomial_terms, coefficients);

		end_time = std :: chrono :: steady_clock::now();

		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		total_time_spent_in_regression += time_span.count();

		vector< vector< vector< datatype > > > weights;
		vector< vector< datatype > > biases;
		system_network.return_network_information(weights, biases);
		vector< vector< vector< vector< datatype > > > > region_descriptions;
		vector< vector < vector< datatype > > > linear_mapping;

		vector< PolynomialApproximator > decomposed_pwls;
		vector< double > lower_bounds;
		vector< double > upper_bounds;

		double tolerance = 1e-4;

		start_time = std :: chrono :: steady_clock::now();

		create_PWL_approximation(polynomial, input_interval, tolerance,
					 region_descriptions, linear_mapping,
			                 decomposed_pwls, lower_bounds,
			                 upper_bounds, linear_piece_count);

		if(linear_piece_count > max_linear_pieces)
		{
			max_linear_pieces = linear_piece_count;
		}

		end_time = std :: chrono :: steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		total_time_spent_in_PWL_construction += time_span.count();

		start_time = std :: chrono :: steady_clock::now();

		vector< datatype > difference;
		system_network.return_interval_difference_wrt_PWL(input_interval, difference, 1, region_descriptions,
			 linear_mapping, offset, scaling, decomposed_pwls, lower_bounds, upper_bounds);

		end_time = std :: chrono :: steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		total_time_spent_in_Sherlock += time_span.count();

		datatype optima = compute_max_abs_in_a_vector(difference);
		if(optima > max_difference)
		{
			max_difference = optima;
		}

		start_time = std :: chrono :: steady_clock::now();

		Polynomial poly_u;

		for(int i=0; i<monomial_terms.size(); ++i)
		{
			monomial_terms[i].insert(monomial_terms[i].begin(), 0);
			monomial_terms[i].push_back(0);

			Monomial monomial(coefficients[i], *((vector<int> *) &monomial_terms[i]));

			poly_u.monomials.push_back(monomial);
		}

		poly_u.reorder();

		// Taylor model approximating the network output

		TaylorModel tm_u;

		vector<Interval> polyRange_initial_set;
		initial_set.tmvPre.polyRangeNormal(polyRange_initial_set, crs.step_end_exp_table);

		poly_u.insert_normal(tm_u, initial_set.tmvPre, polyRange_initial_set, crs.step_end_exp_table, domainDim, crs.cutoff_threshold);

		Interval intError(-optima, optima);
		tm_u.remainder += intError;

		initial_set.tmvPre.tms[2] = tm_u;  // network output becomes new control input

		printf("Step %d\n", k);

		start_time = std :: chrono :: steady_clock::now();

		bool res = plant.reach_symbolic_remainder(result, fp_last, symb_rem, crs, initial_set, 5, 2);

		if(res)
		{
			initial_set = fp_last;
		}
		else
		{
			printf("Terminated due to too large overestimation.\n");
			break;
		}

		end_time = std :: chrono :: steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		total_time_spent_in_Flowstar += time_span.count();
	}

	cout << "Max error encountered = " << max_difference << endl;
	cout << "Max Linear pieces = " << max_linear_pieces << endl;

	std :: chrono :: steady_clock::time_point t2 = std::chrono::steady_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	cout << "Total execution time =  " << time_span.count() << " seconds." << endl;

	cout << "\% of time spent in regression = " << (total_time_spent_in_regression/time_span.count()) * 100.0 << " \% " << endl;
	cout << "\% of  time spent in PWL construction = " << (total_time_spent_in_PWL_construction/time_span.count()) * 100.0<< " \% "  << endl;
	cout << "\% of  time spent in Sherlock = " << (total_time_spent_in_Sherlock/time_span.count()) * 100.0 << " \% " <<  endl;
	cout << "\% of  time spent in Flowstar = " << (total_time_spent_in_Flowstar/time_span.count()) * 100.0 << " \% " << endl;


	// plot the flowpipes in the x-y plane
	FILE *fp = fopen("Plots/SinglePendulum.m", "w");
	plot_2D_interval_MATLAB(fp, "x", "y", stateVars, result);
	fclose(fp);

	save_final_results_to_file(false,4,4,0.01,2,max_difference,time_span.count(),(total_time_spent_in_regression/time_span.count()) * 100.0
	,(total_time_spent_in_PWL_construction/time_span.count()) * 100.0,(total_time_spent_in_Sherlock/time_span.count()) * 100.0,
	(total_time_spent_in_Flowstar/time_span.count()) * 100.0,max_linear_pieces);

	return 0;
}
