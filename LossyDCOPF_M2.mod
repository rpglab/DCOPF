# Lossy DCOPF (Iterative Loss Approximation), 
# developed by Xingpeng.Li
#    Website: https://rpglab.github.io/
# Each iteration: calc total losses and then the change in loss than previous iteration, and 
#    then assign the losses to all/selected buses with same/different participation factors.
# This method converges when the difference of total losses between two consecutive iterations is less than a threshold.
# Possibly, two iterations repeat themselves, which cause method to diverge.

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
param bus_PF	 {BUS}; # Bus Participation Factor, summation over all buses must be 1.

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
###########################


###########################
#### VARIABLES:
var bus_angle  {k in BUS};        # Variable for Bus Angles
var gen_supply {i in GEN};      # Variable for GEN Supply
var line_flow  {j in BRANCH};     # Variable for all line flows
var line_Ploss  {j in BRANCH};     # power loss in each line
var Tloss0;   # system total loss of i-th iteration
var dTloss;   # difference of system total loss between two consecutive periods
var iter;
var PlossBus   {k in BUS};        #
var PlossLine   {j in BRANCH}; 
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

###########################



###########################
#### Load data:
data;
param: BUS: bus_num bus_Pd bus_PF :=       # data for the Buses
        include BusData_LossyDCOPF_M2.dat;

param: GEN: gen_num gen_bus gen_Pmax gen_Pmin gen_C :=     # data for the GENs
        include GenData.dat;

param: BRANCH: branch_num branch_fbus branch_tbus         # data for the Lines 
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
###########################



###########################
option solver gurobi;
option gurobi_options 'mipgap = 0.0';

# option solver cplexamp;
# option cplex_options('mipgap=0.0 integrality = 0.0');

# option solver knitro;
# option knitro_options "outlev=3 alg=1";

let Tloss0 := 0;   # system total loss of i-th iteration
let dTloss := 1;   # difference of system total loss between two consecutive periods
let iter := 0;
repeat while (abs(dTloss) > 1e-9) {
    solve;
    for{j in BRANCH}
    {
        let line_Ploss[j] := 2*branch_Gk[j]*(1-cos(bus_angle[branch_tbus[j]] - bus_angle[branch_fbus[j]]));
	}
	let dTloss := sum{j in BRANCH}line_Ploss[j] - Tloss0;
	let Tloss0 := sum{j in BRANCH}line_Ploss[j];
    for{i in BUS}
    {
		let bus_Pd[i] := bus_Pd[i] + bus_PF[i]*dTloss;
	}	
	let iter := iter+1;
	display iter;
	display dTloss;
};

for{i in BUS}
{
    let PlossBus[i] := 0;
}
for{j in BRANCH}
{
    let PlossLine[j] := 2*branch_r[j]/(branch_r[j]^2+branch_x[j]^2)*(1-cos(bus_angle[branch_fbus[j]] - bus_angle[branch_tbus[j]]));
    let PlossBus[branch_fbus[j]] := if (bus_angle[branch_fbus[j]] > bus_angle[branch_tbus[j]] ) then (PlossBus[branch_fbus[j]] + 0 )
	                                 else if (bus_angle[branch_fbus[j]] <= bus_angle[branch_tbus[j]] ) then (PlossBus[branch_fbus[j]] + PlossLine[j]);
    let PlossBus[branch_tbus[j]] := if (bus_angle[branch_fbus[j]] <= bus_angle[branch_tbus[j]] ) then (PlossBus[branch_tbus[j]] + 0)
	                                 else if (bus_angle[branch_fbus[j]] > bus_angle[branch_tbus[j]] ) then (PlossBus[branch_tbus[j]] + PlossLine[j]);
};

for{i in BRANCH}
{
    let line_PlossTest[i] := 2;
    let line_PlossTest[i] := if (line_flow[i] == branch_rateA[i]) then (1) 
                             else if (line_flow[i] != branch_rateA[i]) then (0);
};

####### Show me the results.  Happy End!

display bus_angle;
display gen_supply;
display line_flow;
display Cost;
display dTloss;
display iter;
display _solve_system_time;
display _solve_user_time;
display _total_solve_time;
option show_stats 1;
display PowerBal.dual;
display PlossBus;
display sum{i in BRANCH}PlossLine[i];
display sum{i in BUS}PlossBus[i];
display line_PlossTest;


