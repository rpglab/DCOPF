# Lossy DCOPF (Piecewise Approximation), developed by Xingpeng.Li
# implemented the model (2)-(14) of the following paper: 
#     O. W. Akinbode and K. W. Hedman "Fictitious losses in the DCOPF with a piecewise linear approximation of losses," 
#     IEEE PES General Meeting, Jul. 2013.


reset;


###########################
#### set
set BUS;    # set of buses
set BRANCH; # set of branches
set GEN;   # Gen Data
###########################


###########################
#### PARAMETERS:
# Bus Data
param bus_num        {BUS}; # Bus Number
param bus_Pd	     {BUS}; # Real Power Demand 

# GEN Data
param gen_num	  {GEN}; # GEN number
param gen_bus	  {GEN}; # GEN location
param gen_Pmax    {GEN}; # Max gen production
param gen_Pmin    {GEN}; # Min gen production when committed
param gen_C	      {GEN}; # Linear Cost Term

# Branch Data
param branch_num     {BRANCH}; # Branch Number
param branch_fbus    {BRANCH}; # from bus for line
param branch_tbus    {BRANCH}; # to bus for line
param branch_r       {BRANCH}; # line resistance
param branch_x       {BRANCH}; # line reactance
param branch_b       {BRANCH}; # total line charging susceptance
param branch_rateA   {BRANCH}; # long term thermal rating
param branch_Gk      {BRANCH}; # Gk = R/(R^2+X^2);

# base MVA
param BaseMVA = 100;  # base MVA
param NumSec = 7;     # number of sections of Piecewise Linear Approximation
param alphas  {k in 1..NumSec}; # alphas of each segment slope for Piecewise Linear Approximation of Power loss
param angleDifMax  {k in 1..NumSec}; #
param angleDifMin  {k in 1..NumSec}; #

###########################


###########################
#### VARIABLES:
var bus_angle  {k in BUS};        # Variable for Bus Angles
var gen_supply {i in GEN};      # Variable for GEN Supply
var line_flow  {j in BRANCH};     # Variable for all line flows
var bus_Ploss  {j in BUS};     # power loss in one line
var line_angleDif  {j in BRANCH};     # angle difference in one line
var line_angleDifP  {j in BRANCH,k in 1..NumSec} >= angleDifMin[k], <=angleDifMax[k];
var line_angleDifN  {j in BRANCH,k in 1..NumSec} >= angleDifMin[k], <=angleDifMax[k];
var line_PlossTest  {j in BRANCH};     # test whether the line is congestion
###########################


###########################
#### OBJECTIVE:
minimize Cost: BaseMVA * sum{i in GEN}gen_supply[i]*gen_C[i];
###########################


###########################
#### CONSTRAINTS:
##--- POWER BALANCE EQUATION.
subject to PowerBal{k in BUS}: 		# Power Balance Constraint
       sum{j in BRANCH: branch_tbus[j] == k}line_flow[j]
	 - sum{j in BRANCH: branch_fbus[j] == k}line_flow[j]
	 - bus_Ploss[k]
	 + sum{i in GEN: gen_bus[i]==k}(gen_supply[i])
	 - bus_Pd[k]   = 0;

	 
##--- POWER OUTPUT LIMITATION FOR GENERATORS.
subject to PGenMaxMin {i in GEN}: gen_Pmin[i] <= gen_supply[i] <= gen_Pmax[i]; # Gen min & max constraint.

##--- LINE FLOW EQUATION, AND THERMAL CONSTRAINTS FOR THE EXISTING LINES.
subject to Line_Flow_value{j in BRANCH}: 
                line_flow[j] = (bus_angle[branch_fbus[j]] - bus_angle[branch_tbus[j]])/branch_x[j];				
subject to Thermal1{j in BRANCH}: 	# Thermal Constraint
	            (-branch_rateA[j]) <= line_flow[j];				
subject to Thermal2{j in BRANCH}:	# Thermal Constraint 2  
	            (branch_rateA[j]) >= line_flow[j];

##--- Power Loss 
subject to Ploss1{j in BRANCH}: 
                (bus_angle[branch_tbus[j]] - bus_angle[branch_fbus[j]])
				   = sum{k in 1..NumSec}line_angleDifP[j,k] - sum{k in 1..NumSec}line_angleDifN[j,k];
subject to Ploss2{k in BUS}: bus_Ploss[k] = 
                sum{j in BRANCH: branch_tbus[j] == k}(2*branch_Gk[j]*(sum{s in 1..NumSec}alphas[s]*line_angleDifN[j,s]))
			  + sum{j in BRANCH: branch_fbus[j] == k}(2*branch_Gk[j]*(sum{s in 1..NumSec}alphas[s]*line_angleDifP[j,s]));
###########################



###########################
#### Load data:
data;
param: BUS: bus_num bus_Pd :=       # data for the Buses
        include BusData.dat;

param: GEN: gen_num gen_bus gen_Pmax gen_Pmin gen_C :=     # data for the GENs
        include GenData.dat;

param: BRANCH: branch_num branch_fbus branch_tbus          # data for the Lines 
               branch_r branch_x branch_b branch_rateA:=
        include BranchData.dat;
###########################



###########################
#----SCALE AND INITIALIZE THE DATA.
for{i in BUS}
{
    let bus_Pd[i] := bus_Pd[i]/BaseMVA;
};
for{i in GEN}
{
    let gen_Pmax[i] := gen_Pmax[i]/BaseMVA;
    let gen_Pmin[i] := gen_Pmin[i]/BaseMVA;
};
for{i in BRANCH}
{
    let branch_rateA[i] := branch_rateA[i]/BaseMVA;
	let branch_Gk[i] := branch_r[i]/(branch_r[i]^2 + branch_x[i]^2);
};
for {k in 1..NumSec}
{
    let alphas[k] := ((1-cos(0.1*k))-(1-cos(0.1*(k-1))))/0.1;
    let angleDifMax[k] := 0.1;
    let angleDifMin[k] := 0;
};
###########################



###########################
option solver gurobi;
option gurobi_options 'mipgap = 0.0';

# option solver cplexamp;
# option cplex_options('mipgap=0.0 integrality = 0.0');

# option solver knitro;
# option knitro_options "outlev=3 alg=1";


solve;

for{i in BRANCH}
{
    let line_PlossTest[i] := 2;
    let line_PlossTest[i] := if (line_flow[i] == branch_rateA[i]) then (1) 
                             else if (line_flow[i] <= branch_rateA[i]) then (0);
};

####### Show me the results.  Happy End!
display _solve_system_time;
display _solve_user_time;
display _total_solve_time;
option show_stats 1;

display Cost;
display bus_angle;
display gen_supply;
display line_flow;
display line_PlossTest;
display PowerBal.dual;  # LMP in $ per 100 MWh energy consumption
display bus_Ploss;
display sum{i in BUS}bus_Ploss[i];
display line_PlossTest;
