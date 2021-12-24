CXX = g++

include Makefile.locale
GUROBI_INCLUDEDIR=$(strip $(GUROBI_PATH))/$(strip $(HOST_ARCH))/include/
GUROBI_LIBDIR=$(strip $(GUROBI_PATH))/$(strip $(HOST_ARCH))/lib/

LIBS = -lgurobi_c++ -lgurobi80 -lflowstar -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk -lmpfi -D_GLIBCXX_USE_CXX11_ABI=0 -m64 -w
CFLAGS = -I . -I ./headers -I /usr/local/include/ -I $(GUROBI_INCLUDEDIR) \
	-I ./headers/eigen_file/ -I ./headers/neural_rule_analysis/src_new -g -O3 -std=c++11
LINK_FLAGS = -g -L ./ -L /usr/local/lib/ -L $(GUROBI_LIBDIR) -L ../flowstar-release -no-pie
OBJS = ./headers/propagate_intervals.o ./headers/network_computation.o ./headers/gurobi_interface.o \
 ./headers/neural_rule_analysis/src_new/AffineArithmeticExpression.o ./headers/neural_rule_analysis/src_new/AffineArithmeticNoiseSymbols.o \
 ./headers/neural_rule_analysis/src_new/Box.o ./headers/neural_rule_analysis/src_new/Monomial.o \
 ./headers/neural_rule_analysis/src_new/mpfiWrapper.o ./headers/neural_rule_analysis/src_new/neuralRuleAnalysisInterfaceMain.o \
 ./headers/neural_rule_analysis/src_new/Polynomial.o ./headers/neural_rule_analysis/src_new/PolynomialApproximator.o \
 ./headers/neural_rule_analysis/src_new/Tiling.o ./headers/neural_rule_analysis/src_new/RangeToVariables.o \
 ./headers/compute_flowpipes.o

all: Ex_1 Ex_2 Ex_3 Ex_4 Ex_5 Ex_6 Ex_7 Ex_8 Ex_9 Ex_10 Ex_Quadrotor Ex_1_LR Ex_2_LR Ex_3_LR Ex_4_LR \
    Ex_5_LR Ex_6_LR Ex_7_LR Ex_8_LR Ex_9_LR Ex_Quadrotor_LR  arch_1 \
  	Ex_1_INT Ex_2_INT Ex_3_INT Ex_4_INT Ex_5_INT Ex_6_INT Ex_7_INT Ex_8_INT Ex_9_INT Ex_Quadrotor_INT \
  	Ex_1_PR Ex_2_PR Ex_3_PR Ex_4_PR Ex_5_PR Ex_6_PR Ex_7_PR Ex_8_PR Ex_9_PR Ex_Quadrotor_PR alter_controller \
  	Ex_10_LR Ex_10_PR Ex_Drone_Model_LR

alter_controller: ./alter_controller.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

Ex_1: ./Ex_1.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_2: ./Ex_2.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_3: ./Ex_3.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_4: ./Ex_4.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_5: ./Ex_5.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_6: ./Ex_6.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_7: ./Ex_7.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_8: ./Ex_8.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_9: ./Ex_9.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_10: ./Ex_10.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_Quadrotor: ./Ex_Quadrotor.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

arch_1: ./arch_1.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)


Ex_1_LR: ./Ex_1_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_2_LR: ./Ex_2_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_3_LR: ./Ex_3_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_4_LR: ./Ex_4_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_5_LR: ./Ex_5_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_6_LR: ./Ex_6_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_7_LR: ./Ex_7_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_8_LR: ./Ex_8_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_9_LR: ./Ex_9_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_Quadrotor_LR: ./Ex_Quadrotor_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

Ex_10_LR: ./Ex_10_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_Drone_Model_LR: ./Ex_Drone_Model_LR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

Ex_10_PR: ./Ex_10_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

Ex_1_INT: ./Ex_1_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_2_INT: ./Ex_2_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_3_INT: ./Ex_3_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_4_INT: ./Ex_4_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_5_INT: ./Ex_5_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_6_INT: ./Ex_6_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_7_INT: ./Ex_7_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_8_INT: ./Ex_8_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_9_INT: ./Ex_9_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_Quadrotor_INT: ./Ex_Quadrotor_INT.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_1_PR: ./Ex_1_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_2_PR: ./Ex_2_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_3_PR: ./Ex_3_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_4_PR: ./Ex_4_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_5_PR: ./Ex_5_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_6_PR: ./Ex_6_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_7_PR: ./Ex_7_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_8_PR: ./Ex_8_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_9_PR: ./Ex_9_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
Ex_Quadrotor_PR: ./Ex_Quadrotor_PR.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

B_Unicycle: ./B_Unicycle.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
B_TORA: ./B_TORA.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
B_ACC: ./B_ACC.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
B_SinglePendulum: ./B_SinglePendulum.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
B_DoublePendulum: ./B_DoublePendulum.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
B_Airplane: ./B_Airplane.o $(OBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)
benchmarks: B_Unicycle B_TORA B_ACC B_SinglePendulum B_DoublePendulum B_Airplane

%.o: %.cc
	$(CXX) -O3 -c $(CFLAGS) -o $@ $< $(LIBS)
%.o: %.cpp
	$(CXX) -O3 -c $(CFLAGS) -o $@ $< $(LIBS)
%.o: %.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f ./headers/*.o *.o Ex_1_LR Ex_2_LR Ex_3_LR Ex_4_LR Ex_5_LR Ex_6_LR Ex_7_LR Ex_8_LR Ex_9_LR Ex_Quadrotor_LR \
		Ex_1_INT Ex_2_INT Ex_3_INT Ex_4_INT Ex_5_INT Ex_6_INT Ex_7_INT Ex_8_INT Ex_9_INT Ex_Quadrotor_INT \
		Ex_1_PR Ex_2_PR Ex_3_PR Ex_4_PR Ex_5_PR Ex_6_PR Ex_7_PR Ex_8_PR Ex_9_PR Ex_Quadrotor_PR  arch_1 \
		Ex_10_LR Ex_10_PR Ex_Drone_Model_LR Ex_1  Ex_2  Ex_3  Ex_4  Ex_5 Ex_6 Ex_7  ./headers/neural_rule_analysis/src_new/*.o
