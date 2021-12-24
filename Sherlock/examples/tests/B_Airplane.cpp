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
	stateVars.declareVar("z");
	stateVars.declareVar("u");
	stateVars.declareVar("v");
	stateVars.declareVar("w");
	stateVars.declareVar("phi");
	stateVars.declareVar("theta");
	stateVars.declareVar("psi");
	stateVars.declareVar("r");
	stateVars.declareVar("p");
	stateVars.declareVar("q");
	stateVars.declareVar("i1");
	stateVars.declareVar("i2");
	stateVars.declareVar("i3");
	stateVars.declareVar("i4");
	stateVars.declareVar("i5");
	stateVars.declareVar("i6");

	int domainDim = 19;

	Expression_AST deriv_x("cos(theta)*cos(psi) * u + (sin(theta)*cos(psi)*sin(phi) - (sin(psi)*cos(phi))) * v + (sin(psi)*sin(phi) + sin(theta)*cos(psi)*cos(phi)) * w", stateVars);
	Expression_AST deriv_y("cos(theta)*sin(psi) * u + (cos(psi)*cos(phi) + sin(theta)*sin(psi)*sin(phi)) * v + (sin(theta)*sin(psi)*cos(phi) - (cos(psi)*sin(phi))) * w", stateVars);
	Expression_AST deriv_z("-sin(theta) * u + cos(theta)*sin(phi) * v + cos(theta)*cos(phi) * w", stateVars);
	Expression_AST deriv_u("-sin(theta) + i1 - q * w + r * v", stateVars);
	Expression_AST deriv_v("cos(theta) * sin(phi) + i2 - r * u + p * w", stateVars);
	Expression_AST deriv_w("cos(theta) * cos(phi) + i3 - p * v + q * u", stateVars);
	Expression_AST deriv_phi("p + sin(theta)/cos(theta) * sin(phi) * q + sin(theta)/cos(theta) * cos(phi) * r", stateVars);
	Expression_AST deriv_theta("cos(phi) * q - sin(phi) * r", stateVars);
	Expression_AST deriv_psi("1 / cos(theta) * sin(phi) * q + 1 / cos(theta) * cos(phi) * r", stateVars);
	Expression_AST deriv_r("i6", stateVars);
	Expression_AST deriv_p("i5", stateVars);
	Expression_AST deriv_q("i4", stateVars);
	Expression_AST deriv_i1("0", stateVars);
	Expression_AST deriv_i2("0", stateVars);
	Expression_AST deriv_i3("0", stateVars);
	Expression_AST deriv_i4("0", stateVars);
	Expression_AST deriv_i5("0", stateVars);
	Expression_AST deriv_i6("0", stateVars);

	ODE plant(stateVars);
	plant.assignDerivative("x", deriv_x);
	plant.assignDerivative("y", deriv_y);
	plant.assignDerivative("z", deriv_z);
	plant.assignDerivative("u", deriv_u);
	plant.assignDerivative("v", deriv_v);
	plant.assignDerivative("w", deriv_w);
	plant.assignDerivative("phi", deriv_phi);
	plant.assignDerivative("theta", deriv_theta);
	plant.assignDerivative("psi", deriv_psi);
	plant.assignDerivative("r", deriv_r);
	plant.assignDerivative("p", deriv_p);
	plant.assignDerivative("q", deriv_q);
	plant.assignDerivative("i1", deriv_i1);
	plant.assignDerivative("i2", deriv_i2);
	plant.assignDerivative("i3", deriv_i3);
	plant.assignDerivative("i4", deriv_i4);
	plant.assignDerivative("i5", deriv_i5);
	plant.assignDerivative("i6", deriv_i6);


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
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	estimation.push_back(E);
	crs.setRemainderEstimation(estimation);

	// call this function whenever a parameter is set or changed
	crs.prepareForReachability();

	// Simple range propagation
	char controller_file[] = "controller_Airplane_sherlock";

	network_handler system_network(controller_file);

	vector< vector< vector< datatype > > > buff_weights;
	vector< vector < datatype > > buff_biases;
	system_network.return_network_information(buff_weights, buff_biases);
	system_network.cast_to_single_output_network(buff_weights, buff_biases, 1);
	network_handler control_output_1(buff_weights, buff_biases);

	system_network.return_network_information(buff_weights, buff_biases);
	system_network.cast_to_single_output_network(buff_weights, buff_biases, 2);
	network_handler control_output_2(buff_weights, buff_biases);

	system_network.return_network_information(buff_weights, buff_biases);
	system_network.cast_to_single_output_network(buff_weights, buff_biases, 3);
	network_handler control_output_3(buff_weights, buff_biases);

	system_network.return_network_information(buff_weights, buff_biases);
	system_network.cast_to_single_output_network(buff_weights, buff_biases, 4);
	network_handler control_output_4(buff_weights, buff_biases);

	system_network.return_network_information(buff_weights, buff_biases);
	system_network.cast_to_single_output_network(buff_weights, buff_biases, 5);
	network_handler control_output_5(buff_weights, buff_biases);

	system_network.return_network_information(buff_weights, buff_biases);
	system_network.cast_to_single_output_network(buff_weights, buff_biases, 6);
	network_handler control_output_6(buff_weights, buff_biases);

	network_handler control_outputs[] = {control_output_1, control_output_2,
		control_output_3, control_output_4, control_output_5, control_output_6
	};

	/*
	 * Initial set can be a box which is represented by a vector of intervals.
	 * The i-th component denotes the initial set of the i-th state variable.
	 */
// 	Interval init_x(0.0,0.00001), init_y(0.0,0.00001), init_z(0.0,0.00001),  // full set
// 		 init_u(0.0,1.0), init_v(0.0,1.0), init_w(0.0,1.0),
// 		 init_phi(0.0,1.0), init_theta(0.0,1.0), init_psi(0.0,1.0),
// 		 init_r(0.0,0.00001), init_p(0.0,0.00001), init_q(0.0,0.00001),
// 	         init_i1, init_i2, init_i3, init_i4, init_i5, init_i6, intZero;
// 	Interval init_x(0.0,0.01), init_y(0.0,0.01), init_z(0.0,0.01),  // small subset
// 		 init_u(0.99,1.0), init_v(0.99,1.0), init_w(0.99,1.0),
// 		 init_phi(0.99,1.0), init_theta(0.99,1.0), init_psi(0.99,1.0),
// 		 init_r(0.0,0.01), init_p(0.0,0.01), init_q(0.0,0.01),
	Interval init_x(0.0,0.00001), init_y(0.0,0.00001), init_z(0.0,0.00001),  // very small subset
		 init_u(0.99999,1.0), init_v(0.99999,1.0), init_w(0.99999,1.0),
		 init_phi(0.99999,1.0), init_theta(0.99999,1.0), init_psi(0.99999,1.0),
		 init_r(0.0,0.00001), init_p(0.0,0.00001), init_q(0.0,0.00001),
	         init_i1, init_i2, init_i3, init_i4, init_i5, init_i6, intZero;
	std::vector<Interval> X0;
	X0.push_back(init_x);
	X0.push_back(init_y);
	X0.push_back(init_z);
	X0.push_back(init_u);
	X0.push_back(init_v);
	X0.push_back(init_w);
	X0.push_back(init_phi);
	X0.push_back(init_theta);
	X0.push_back(init_psi);
	X0.push_back(init_r);
	X0.push_back(init_p);
	X0.push_back(init_q);
	X0.push_back(init_i1);
	X0.push_back(init_i2);
	X0.push_back(init_i3);
	X0.push_back(init_i4);
	X0.push_back(init_i5);
	X0.push_back(init_i6);


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

		vector< vector< datatype > > input_interval(12, vector< datatype >(2,0));
		for (int i = 0; i < 12; ++i) {
			input_interval[i][0] = NN_input[i].inf();
			input_interval[i][1] = NN_input[i].sup();
		}

		printf("[%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf], \t [%lf, %lf]\n",
		       input_interval[0][0], input_interval[0][1],
		       input_interval[1][0], input_interval[1][1],
		       input_interval[2][0], input_interval[2][1],
		       input_interval[3][0], input_interval[3][1],
		       input_interval[4][0], input_interval[4][1],
		       input_interval[5][0], input_interval[5][1],
		       input_interval[6][0], input_interval[6][1],
		       input_interval[7][0], input_interval[7][1],
		       input_interval[8][0], input_interval[8][1],
		       input_interval[9][0], input_interval[9][1],
		       input_interval[10][0], input_interval[10][1],
		       input_interval[11][0], input_interval[11][1]);

		vector< vector< datatype > > input_constraints;
		create_constraint_from_interval(input_constraints, input_interval);

		datatype offset = 0.0;
		datatype scaling = 1.0;
		unsigned int degree = 2;
		
		for (int o = 0; o < 6; ++o) {
			network_handler control_output = control_outputs[o];

			vector< vector< unsigned int > > monomial_terms;
			vector< datatype > coefficients;
			start_time = std :: chrono :: steady_clock::now();
			generate_polynomial_for_NN(control_output, degree, input_constraints,
						offset, scaling, monomial_terms, coefficients);
			int i,j,l;
			vector<my_monomial_t> polynomial;
			polynomial = create_polynomial_from_monomials_and_coeffs(monomial_terms, coefficients);

			end_time = std :: chrono :: steady_clock::now();

			time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
			total_time_spent_in_regression += time_span.count();

			vector< vector< vector< datatype > > > weights;
			vector< vector< datatype > > biases;
			control_output.return_network_information(weights, biases);
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
			control_output.return_interval_difference_wrt_PWL(input_interval, difference, 1,
				region_descriptions, linear_mapping, offset, scaling, decomposed_pwls,
				lower_bounds, upper_bounds);

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

			poly_u.insert_normal(tm_u, initial_set.tmvPre, polyRange_initial_set,
					     crs.step_end_exp_table, domainDim, crs.cutoff_threshold);

			Interval intError(-optima, optima);
			tm_u.remainder += intError;

			initial_set.tmvPre.tms[12+o] = tm_u;  // network output becomes new control input
		}

		printf("Step %d\n", k);

		start_time = std :: chrono :: steady_clock::now();

		bool res = plant.reach_symbolic_remainder(result, fp_last, symb_rem, crs, initial_set, 10, 2);

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


	// plot the flowpipes
	FILE *fp = fopen("Plots/Airplane_yphi.m", "w");
	plot_2D_interval_MATLAB(fp, "y", "phi", stateVars, result);
	fclose(fp);
	FILE *fp2 = fopen("Plots/Airplane_thetapsi.m", "w");
	plot_2D_interval_MATLAB(fp2, "theta", "psi", stateVars, result);
	fclose(fp2);

	save_final_results_to_file(false,6,4,0.01,2,max_difference,time_span.count(),(total_time_spent_in_regression/time_span.count()) * 100.0
	,(total_time_spent_in_PWL_construction/time_span.count()) * 100.0,(total_time_spent_in_Sherlock/time_span.count()) * 100.0,
	(total_time_spent_in_Flowstar/time_span.count()) * 100.0,max_linear_pieces);

	return 0;
}
