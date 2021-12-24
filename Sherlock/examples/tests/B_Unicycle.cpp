// original file: Ex_10_PR

#include "./headers/propagate_intervals.h"
#include "../flowstar-release/Continuous.h"

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

	/*
	 * Declaration of the state variables.
	 * The first one should always be the local time variable which can be
	 * viewed as a preserved variable. It is only used internally by the library.
	 */
	stateVars.declareVar("t");
	stateVars.declareVar("x1");
	stateVars.declareVar("x2");
  stateVars.declareVar("x3");
  stateVars.declareVar("x4");
	stateVars.declareVar("u1");
  stateVars.declareVar("u2");

	int domainDim = 7;

	// Expression_AST deriv_x1( " x4*cos(x3 + 0.4545*(sin(u2)/cos(u2)) - ((0.4545*(sin(u2)/cos(u2)))^3)/3 + ((0.4545*(sin(u2)/cos(u2)))^5)/5) " , stateVars);
  // Expression_AST deriv_x2( " x4*sin(x3 + 0.4545*(sin(u2)/cos(u2)) - ((0.4545*(sin(u2)/cos(u2)))^3)/3 + ((0.4545*(sin(u2)/cos(u2)))^5)/5) " , stateVars);
  // Expression_AST deriv_x3( " (x4 * sin( 0.4545*(sin(u2)/cos(u2)) - ((0.4545*(sin(u2)/cos(u2)))^3)/3 + ((0.4545*(sin(u2)/cos(u2)))^5)/5) )/1.5 ", stateVars);

  Expression_AST deriv_x1( " x4 * cos(x3) " , stateVars);
  Expression_AST deriv_x2( " x4 * sin(x3) " , stateVars);
  Expression_AST deriv_x3( "u2 ", stateVars);
  Expression_AST deriv_x4("u1 + [-0.0001,0.0001] ", stateVars);
  Expression_AST deriv_u1("0", stateVars);
  Expression_AST deriv_u2("0", stateVars);



	ODE plant(stateVars);
	plant.assignDerivative("x1", deriv_x1);
	plant.assignDerivative("x2", deriv_x2);
  plant.assignDerivative("x3", deriv_x3);
	plant.assignDerivative("x4", deriv_x4);
	plant.assignDerivative("u1", deriv_u1);
  plant.assignDerivative("u2", deriv_u2);


	/*
	 * Specify the parameters for reachability computation.
	 */
	Continuous_Reachability_Setting crs;

	// step size
	crs.setFixedStepsize(0.01);

	// Taylor model order
	crs.setFixedOrder(30);

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
  estimation.push_back(E);	// estimation for the 6th variable
  estimation.push_back(E);	// estimation for the 7th variable
	crs.setRemainderEstimation(estimation);

	// call this function whenever a parameter is set or changed
	crs.prepareForReachability();





	// Simple range propagation
// 	char controller_file[] = "../systems_with_networks/Ex_car_model/NN_sherlock/neural_network_controller_2" ;
	char controller_file[] = "controller_Unicycle_sherlock";

	network_handler system_network(controller_file);

  vector< vector< vector< datatype > > > buff_weights;
  vector< vector < datatype > > buff_biases;
  system_network.return_network_information(buff_weights, buff_biases);
  system_network.cast_to_single_output_network(buff_weights, buff_biases, 1);
  network_handler control_output_1(buff_weights, buff_biases);

  system_network.return_network_information(buff_weights, buff_biases);
  system_network.cast_to_single_output_network(buff_weights, buff_biases, 2);
  network_handler control_output_2(buff_weights, buff_biases);


	/*
	 * Initial set can be a box which is represented by a vector of intervals.
	 * The i-th component denotes the initial set of the i-th state variable.
	 */
	 Interval init_x1(9.5,9.55), init_x2(-4.5,-4.45), init_x3(2.1,2.11), init_x4(1.5,1.51), init_u1, init_u2, intZero;
	std::vector<Interval> X0;
	X0.push_back(init_x1);
	X0.push_back(init_x2);
  X0.push_back(init_x3);
	X0.push_back(init_x4);
	X0.push_back(init_u1);
  X0.push_back(init_u2);


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
	int no_of_piecewise_linear_pieces= 0 ;
	int linear_piece_count = 0;
	int max_linear_pieces = 0;

	for(int k=0; k<50; ++k)
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


		printf("[%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf]\n", input_interval[0][0], input_interval[0][1],
				input_interval[1][0], input_interval[1][1], input_interval[2][0], input_interval[2][1], input_interval[3][0], input_interval[3][1]);

		vector< vector< datatype > > input_constraints;
		create_constraint_from_interval(input_constraints, input_interval);


    // /////////  Control input 1 /////////

		vector< vector< unsigned int > > monomial_terms_1;
		vector< datatype > coefficients_1;
		datatype offset = -20.0;
    // datatype offset = -100.0;

		datatype scaling = 1.0;
		unsigned int degree = 2;

		start_time = std :: chrono :: steady_clock::now();

		generate_polynomial_for_NN(control_output_1, degree, input_constraints, offset, scaling, monomial_terms_1, coefficients_1);
		int i,j,l;

		vector<my_monomial_t> polynomial_1;
		polynomial_1 = create_polynomial_from_monomials_and_coeffs(monomial_terms_1, coefficients_1);

		end_time = std :: chrono :: steady_clock::now();

		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		total_time_spent_in_regression += time_span.count();

		vector< vector< vector< datatype > > > weights;
		vector< vector< datatype > > biases;
		control_output_1.return_network_information(weights, biases);

		vector< vector< vector< vector< datatype > > > > region_descriptions_1;
		vector< vector < vector< datatype > > > linear_mapping_1;
    //
    //
		vector< PolynomialApproximator > decomposed_pwls_1;
		vector< double > lower_bounds_1;
		vector< double > upper_bounds_1;

		double tolerance = 1e-5;

		start_time = std :: chrono :: steady_clock::now();

		create_PWL_approximation(polynomial_1, input_interval, tolerance, region_descriptions_1, linear_mapping_1,
			 decomposed_pwls_1, lower_bounds_1, upper_bounds_1, linear_piece_count);


			 if(linear_piece_count > max_linear_pieces)
 	 		{
 	 			max_linear_pieces = linear_piece_count;
 	 		}


	  end_time = std :: chrono :: steady_clock::now();
 	  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
 	  total_time_spent_in_PWL_construction += time_span.count();


	  start_time = std :: chrono :: steady_clock::now();


		vector< datatype > difference;
		control_output_1.return_interval_difference_wrt_PWL(input_interval, difference, 1, region_descriptions_1,
			 linear_mapping_1, offset, scaling, decomposed_pwls_1, lower_bounds_1, upper_bounds_1);

	 end_time = std :: chrono :: steady_clock::now();
  	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		total_time_spent_in_Sherlock += time_span.count();


		datatype optima = compute_max_abs_in_a_vector(difference);
		// cout << " Difference is : " << optima << endl;
		if(optima > max_difference)
		{
			max_difference = optima;
		}

		start_time = std :: chrono :: steady_clock::now();

		Polynomial poly_u_1;

		for(int i=0; i < monomial_terms_1.size(); ++i)
		{
			monomial_terms_1[i].insert(monomial_terms_1[i].begin(), 0);
			monomial_terms_1[i].push_back(0);

			Monomial monomial(coefficients_1[i], *((vector<int> *) & monomial_terms_1[i]));

			poly_u_1.monomials.push_back(monomial);
		}

		poly_u_1.reorder();

		TaylorModel tm_u_1;

		vector<Interval> polyRange_initial_set_1;
		initial_set.tmvPre.polyRangeNormal(polyRange_initial_set_1, crs.step_end_exp_table);

		poly_u_1.insert_normal(tm_u_1, initial_set.tmvPre, polyRange_initial_set_1, crs.step_end_exp_table, domainDim, crs.cutoff_threshold);

		Interval intError(-optima, optima);
		tm_u_1.remainder += intError;

		initial_set.tmvPre.tms[4] = tm_u_1;

    /////// Control input 2 ///////

    vector< vector< unsigned int > > monomial_terms_2;
		vector< datatype > coefficients_2;
		// datatype offset = -100.0;
		// datatype scaling = 1.0;
		// unsigned int degree = 1;

		start_time = std :: chrono :: steady_clock::now();

		generate_polynomial_for_NN(control_output_2, degree, input_constraints, offset, scaling, monomial_terms_2, coefficients_2);
		// int i,j,l;

    // i = 0;
    // for(auto coeff : coefficients_2)
    // {
    //   cout << coeff << " * [ " << monomial_terms_2[i][0] <<monomial_terms_2[i][1] << monomial_terms_2[i][2] << monomial_terms_2[i][3] << " ] " ;
    //   i++;
    // }
    // cout << endl;

		vector<my_monomial_t> polynomial_2;
		polynomial_2 = create_polynomial_from_monomials_and_coeffs(monomial_terms_2, coefficients_2);

		end_time = std :: chrono :: steady_clock::now();

		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		total_time_spent_in_regression += time_span.count();

		// vector< vector< vector< datatype > > > weights;
		// vector< vector< datatype > > biases;
    weights.clear();
    biases.clear();
		control_output_2.return_network_information(weights, biases);

		vector< vector< vector< vector< datatype > > > > region_descriptions_2;
		vector< vector < vector< datatype > > > linear_mapping_2;


		vector< PolynomialApproximator > decomposed_pwls_2;
		vector< double > lower_bounds_2;
		vector< double > upper_bounds_2;

		tolerance = 1e-5;

		start_time = std :: chrono :: steady_clock::now();

		create_PWL_approximation(polynomial_2, input_interval, tolerance, region_descriptions_2, linear_mapping_2,
			 decomposed_pwls_2, lower_bounds_2, upper_bounds_2, linear_piece_count);

			 if(linear_piece_count > max_linear_pieces)
		 	{
		 		max_linear_pieces = linear_piece_count;
		 	}

		 end_time = std :: chrono :: steady_clock::now();
	 	 time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
	 	 total_time_spent_in_PWL_construction += time_span.count();


		 start_time = std :: chrono :: steady_clock::now();

		difference.clear();
		control_output_2.return_interval_difference_wrt_PWL(input_interval, difference, 1, region_descriptions_2,
			 linear_mapping_2, offset, scaling, decomposed_pwls_2, lower_bounds_2, upper_bounds_2);

		 end_time = std :: chrono :: steady_clock::now();
	  	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
			total_time_spent_in_Sherlock += time_span.count();

		optima = compute_max_abs_in_a_vector(difference);
		// cout << " Difference is : " << optima << endl;

		if(optima > max_difference)
		{
			max_difference = optima;
		}

		Polynomial poly_u_2;

		for(int i=0; i<monomial_terms_2.size(); ++i)
		{
			monomial_terms_2[i].insert(monomial_terms_2[i].begin(), 0);
			monomial_terms_2[i].push_back(0);

			Monomial monomial(coefficients_2[i], *((vector<int> *) &monomial_terms_2[i]));

			poly_u_2.monomials.push_back(monomial);
		}

		poly_u_2.reorder();

		TaylorModel tm_u_2;

		vector<Interval> polyRange_initial_set_2;
		initial_set.tmvPre.polyRangeNormal(polyRange_initial_set_2, crs.step_end_exp_table);

/*
		for(int i=0; i<initial_set.tmvPre.tms.size(); ++i)
		{
			initial_set.tmvPre.tms[i].dump_interval(stdout, stateVars.varNames);printf("\n");
		}
		printf("\n");
*/

		poly_u_2.insert_normal(tm_u_2, initial_set.tmvPre, polyRange_initial_set_2, crs.step_end_exp_table, domainDim, crs.cutoff_threshold);

		Interval intError_2(-optima, optima);
		tm_u_2.remainder += intError_2;

		initial_set.tmvPre.tms[5] = tm_u_2;

    ////////////////////////////////////////

		printf("Step %d\n", k);

		start_time = std :: chrono :: steady_clock::now();



		bool res = plant.reach_symbolic_remainder(result, fp_last, symb_rem, crs, initial_set, 20, 200);
		// bool res = plant.reach_interval_remainder(result, fp_last, crs, initial_set, 10);

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
	FILE *fp = fopen("Plots/Unicycle12.m", "w");
	plot_2D_interval_MATLAB(fp, "x1", "x2", stateVars, result);
	fclose(fp);
	FILE *fp2 = fopen("Plots/Unicycle34.m", "w");
	plot_2D_interval_MATLAB(fp2, "x3", "x4", stateVars, result);
	fclose(fp2);

	save_final_results_to_file(false,1,30,0.01,2,max_difference,time_span.count(),(total_time_spent_in_regression/time_span.count()) * 100.0
	,(total_time_spent_in_PWL_construction/time_span.count()) * 100.0,(total_time_spent_in_Sherlock/time_span.count()) * 100.0,
	(total_time_spent_in_Flowstar/time_span.count()) * 100.0,max_linear_pieces);




  return 0;
}
