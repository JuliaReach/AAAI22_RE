// original file: Ex_10_PR

#include "./headers/propagate_intervals.h"
#include "../flowstar-release/Continuous.h"
#include "neuralRuleAnalysisInterface.h"

using namespace std;
using namespace flowstar;

datatype offset_in_constraint_comb = constr_comb_offset;
// So the above data assumes, that this number is
// same for all the networks you have in the setting,
// and also all the networks are single output network
// with the first one giving the output

int main()
{
	Variables stateVars;
	bool print_polynomial = false;

	/*
	 * Declaration of the state variables.
	 * The first one should always be the local time variable which can be
	 * viewed as a preserved variable. It is only used internally by the library.
	 */
	stateVars.declareVar("t");
	stateVars.declareVar("x");
	stateVars.declareVar("y");
	stateVars.declareVar("z");
	stateVars.declareVar("w");
	stateVars.declareVar("u");

	int domainDim = 6;


	Expression_AST deriv_x("y", stateVars);
	Expression_AST deriv_y("-x + 0.1*sin(z)", stateVars);
	Expression_AST deriv_z("w", stateVars);
	Expression_AST deriv_w("u", stateVars);
	Expression_AST deriv_u("0", stateVars);


	ODE plant(stateVars);
	plant.assignDerivative("x", deriv_x);
	plant.assignDerivative("y", deriv_y);
	plant.assignDerivative("z", deriv_z);
	plant.assignDerivative("w", deriv_w);
	plant.assignDerivative("u", deriv_u);


	/*
	 * Specify the parameters for reachability computation.
	 */
	Continuous_Reachability_Setting crs;

	// step size
	if(print_polynomial)
	{
		crs.setFixedStepsize(0.01);
	}
	else
	{
		crs.setFixedStepsize(0.1);
	}

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
	estimation.push_back(E);	// estimation for the 4th variable
	estimation.push_back(E);	// estimation for the 5th variable
	crs.setRemainderEstimation(estimation);

	// call this function whenever a parameter is set or changed
	crs.prepareForReachability();





	// Simple range propagation
	// char controller_file[] = "../systems_with_networks/Ex_Tora/neural_network_controller_main" ;
// 	char controller_file[] = "../systems_with_networks/Ex_Tora/bigger_controller" ;
	char controller_file[] = "controller_TORA_sherlock";

	network_handler system_network(controller_file);



	/*
	 * Initial set can be a box which is represented by a vector of intervals.
	 * The i-th component denotes the initial set of the i-th state variable.
	 */
	// Interval init_x(0.5,0.51), init_y(-0.4,-0.39), init_z(-0.1,-0.09), init_w(0.1,0.12), init_u, intZero;
	// Interval init_x(0.6,0.9), init_y(-0.7,-0.4), init_z(-0.4,-0.1), init_w(0.5,0.8), init_u, intZero;
	Interval init_x(0.6,0.7), init_y(-0.7,-0.6), init_z(-0.4,-0.3), init_w(0.5,0.6), init_u, intZero;
  // Interval init_x(0.8,0.85), init_y(-0.5,-0.45), init_z(-0.2,-0.15), init_w(0.7,0.75), init_u, intZero;
	std::vector<Interval> X0;
	X0.push_back(init_x);
	X0.push_back(init_y);
	X0.push_back(init_z);
	X0.push_back(init_w);
	X0.push_back(init_u);


	// translate the initial set to a flowpipe
	Flowpipe initial_set(X0, intZero);

	// the flowpipe that keeps the overapproximation at the end of a time horizon
	Flowpipe fp_last;

	// the symbolic remainder
	Symbolic_Remainder symb_rem(initial_set);

	std::list<Flowpipe> result;
	std::list<Flowpipe> flowpipes_end;

	flowpipes_end.push_back(initial_set);
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
	int linear_piece_count = 0;
	int max_linear_pieces = 0;

	for(int k=0; k<20; ++k)
	{
		std::vector<Interval> NN_input;
		initial_set.intEvalNormal(NN_input, crs.step_end_exp_table, crs.cutoff_threshold);

		vector< vector< datatype > > input_interval(4, vector< datatype >(2,0));
		input_interval[0][0] = NN_input[0].inf();
		input_interval[0][1] = NN_input[0].sup();

		input_interval[1][0] = NN_input[1].inf();
		input_interval[1][1] = NN_input[1].sup();

		input_interval[2][0] = NN_input[2].inf();
		input_interval[2][1] = NN_input[2].sup();

		input_interval[3][0] = NN_input[3].inf();
		input_interval[3][1] = NN_input[3].sup();

		printf("[%lf, %lf], \t [%lf, %lf]\t [%lf, %lf]\t [%lf, %lf]\n", input_interval[0][0], input_interval[0][1],
				input_interval[1][0], input_interval[1][1], input_interval[2][0], input_interval[2][1], input_interval[3][0], input_interval[3][1]);


		vector< vector< datatype > > input_constraints;
		create_constraint_from_interval(input_constraints, input_interval);


		vector< vector< unsigned int > > monomial_terms;
		vector< datatype > coefficients;
		datatype offset = -10;
		datatype scaling = 1;


		unsigned int degree = 2;

		start_time = std :: chrono :: steady_clock::now();

		generate_polynomial_for_NN(system_network, degree, input_constraints, offset, scaling, monomial_terms, coefficients);
		int i,j,l;

		if(print_polynomial)
		{
			cout << " The polynomial we are trying to fit is the following : " << endl;
			for(int i = 0; i < monomial_terms.size(); i++)
			{
				cout << "For monomial term [ " ;
				for(int j = 0; j < monomial_terms[i].size(); j++)
				{
					cout << " , " << monomial_terms[i][j] ;
				}
				cout << " ]  ------ " << coefficients[i] << endl;
			}
		}

		vector<my_monomial_t> polynomial;
		polynomial = create_polynomial_from_monomials_and_coeffs(monomial_terms, coefficients);

		end_time = std :: chrono :: steady_clock::now();

		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		total_time_spent_in_regression += time_span.count();

		vector< vector< vector< datatype > > > weights;
		vector< vector< datatype > > biases;
		system_network.return_network_information(weights, biases);

		// vector< vector< vector< vector< datatype > > > > region_descriptions;
		// vector< vector < vector< datatype > > > linear_mapping;
    // //
    // //
		// double tolerance = 1e-5;
		// create_PWL_approximation(polynomial, input_interval, tolerance, region_descriptions, linear_mapping);
		//
		//
		// vector< datatype > difference;
		// system_network.return_interval_difference_wrt_PWL(input_interval, difference, 1, region_descriptions, linear_mapping, offset, scaling);
		//
		// datatype optima = compute_max_abs_in_a_vector(difference);
		// cout << "Difference is : " << optima << endl;

		vector< vector< vector< vector< datatype > > > > region_descriptions;
		vector< vector < vector< datatype > > > linear_mapping;
    //
    //
		vector< PolynomialApproximator > decomposed_pwls;
		vector< double > lower_bounds;
		vector< double > upper_bounds;
//
		double tolerance = 1e-4;

		start_time = std :: chrono :: steady_clock::now();

		create_PWL_approximation(polynomial, input_interval, tolerance, region_descriptions, linear_mapping,
			 decomposed_pwls, lower_bounds, upper_bounds, linear_piece_count);

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
		if(print_polynomial)
		{
			cout << "Max difference between polynomial and the neural network is --- " << optima << endl;
		}
		// cout << " Difference is : " << optima << endl;

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

		TaylorModel tm_u;

		vector<Interval> polyRange_initial_set;
		initial_set.tmvPre.polyRangeNormal(polyRange_initial_set, crs.step_end_exp_table);



		poly_u.insert_normal(tm_u, initial_set.tmvPre, polyRange_initial_set, crs.step_end_exp_table, domainDim, crs.cutoff_threshold);

		Interval intError(-optima, optima);
		tm_u.remainder += intError;

		initial_set.tmvPre.tms[4] = tm_u;

		printf("Step %d\n", k);




		bool res = plant.reach_symbolic_remainder(result, fp_last, symb_rem, crs, initial_set, 10, 200);
//		bool res = plant.reach_interval_remainder(result, fp_last, crs, initial_set, 10);

		if(res)
		{
			initial_set = fp_last;

//			flowpipes_end.push_back(fp_last);
		}
		else
		{
			printf("Terminated due to too large overestimation.\n");
			break;
		}


		end_time = std :: chrono :: steady_clock::now();
  	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		total_time_spent_in_Flowstar += time_span.count();

		if(print_polynomial)
		{
			cout << "Next set computed  : " << endl;
			std::vector<Interval> NN_input;
			fp_last.intEvalNormal(NN_input, crs.step_end_exp_table, crs.cutoff_threshold);

			vector< vector< datatype > > input_interval(4, vector< datatype >(2,0));

			input_interval[0][0] = NN_input[0].inf();
			input_interval[0][1] = NN_input[0].sup();

			input_interval[1][0] = NN_input[1].inf();
			input_interval[1][1] = NN_input[1].sup();

			input_interval[2][0] = NN_input[2].inf();
			input_interval[2][1] = NN_input[2].sup();

			input_interval[3][0] = NN_input[3].inf();
			input_interval[3][1] = NN_input[3].sup();

			printf("[%lf, %lf], \t [%lf, %lf]\t [%lf, %lf]\t [%lf, %lf]\n", input_interval[0][0], input_interval[0][1],
					input_interval[1][0], input_interval[1][1], input_interval[2][0], input_interval[2][1], input_interval[3][0], input_interval[3][1]);

			break;
		}


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
	FILE *fp = fopen("Plots/TORA12.m", "w");
	plot_2D_interval_MATLAB(fp, "x", "y", stateVars, result);
	fclose(fp);
	FILE *fp2 = fopen("Plots/TORA34.m", "w");
	plot_2D_interval_MATLAB(fp2, "z", "w", stateVars, result);
	fclose(fp2);



	save_final_results_to_file(false,2,4,0.1,2,max_difference,time_span.count(),(total_time_spent_in_regression/time_span.count()) * 100.0
	,(total_time_spent_in_PWL_construction/time_span.count()) * 100.0,(total_time_spent_in_Sherlock/time_span.count()) * 100.0,
	(total_time_spent_in_Flowstar/time_span.count()) * 100.0,max_linear_pieces);



	return 0;
}
