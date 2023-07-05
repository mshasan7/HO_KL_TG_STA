using PyCall
@pyimport matplotlib
matplotlib.use("Agg")

using PyPlot
using DelimitedFiles
using Printf 
using FFTW
using HDF5
#--------------load the functionfiles-------------------------------------
include("initialise.jl")
include("functions.jl")
#------------------------------
hbar =1.0
m=1.0;
#----------------------------
mydict = Dict{String,Any}(
"nx" =>2^15,
"omega_0"=>(1.),
"omega_f"=>(8.), 
"N_bec"=>1,
"tolerance"=>1e-13,
"totspace"=>20,
"imaginary_time"=>true,
"real_time"=>true,
"images_output"=>true,
"save_wfc"=>true,
"N_pic"=>100,
);
if(haskey(mydict, "nx")) mydict["dx"] = mydict["totspace"]/mydict["nx"] end;
mydict["g"] = (hbar.^2*pi^2) ./(2.0.*m);
mydict["gamma"] = sqrt((mydict["omega_0"]./mydict["omega_f"]));
#---------------------------------------------------------------------
prt = Param(mydict);
opr = Operators(mydict);
#-----------------------------------
opr.V = potential_w0(prt);
opr.V_wf = potential_wf(prt);
opr.wfc = bec_wavefunction(prt);
opr.wfc_wf = bec_wavefunction(prt);
opr.wfc_in .= opr.wfc;
#----------------varying tf-----------------
tf_vec = collect(0.1:0.1:1)
@show ARGS
my_argument = parse(Int32,ARGS[1])
opr.tf = tf_vec[my_argument];
#----------------------------GS calculation------------------------
if(prt.im_time)
    opr.im_value = 1.0;
    opr.K = exp.(-opr.im_value.*0.25*(prt.kx.^2)*opr.dt);
    opr.K_wf = exp.(-opr.im_value.*0.25*(prt.kx.^2)*opr.dt);
    @time split_method!(prt, opr);
else
    println("imaginary time parameter off")
end
#----------------------save the GS wfcs--------------------------
#if(prt.save_wfc)
#    h5open("data/opr.wfc_gs_omega_0_tf_.h5","w") do h5
#               h5["opr.wfc"]= opr.wfc
#    end;
#    #h5open("data/opr.wfc_gs_omega_f_tf_.h5","w") do h5
#    #           h5["opr.wfc_wf"]= opr.wfc_wf
#    #end;
#else
#    println("save parameter is off")
#end
#------------------------ground state calculation ends-----------------------
gs_plots!(prt,opr);
PyPlot.savefig("data/ground_state_density.png", bbox_inches = "tight", pad_inches = 0.1)
PyPlot.clf();
#----------------define real-time parameters----------------
opr.dt = 5e-7;
#------------------time evolution begins----------------------------------------
if(prt.re_time)
    opr.im_value = 1.0.*im;
    @. opr.K = exp.(-opr.im_value .* 0.25.*(prt.kx.^2).*opr.dt);
    split_method!(prt, opr);
else
    println("real time parameter is switched off")
end
#-------------------save the wfcs-------------------
if(prt.save_wfc)
    h5open("data/opr.wfc_t_tf"*string(opr.tf)*".h5","w") do h5
               h5["opr.wfc_t"]= opr.wfc_t;
    end
else
    println("save parameter is off")
end
rt_plots!(prt,opr)
PyPlot.savefig("data/rte_sta_"*string(opr.tf)*"_.png", bbox_inches = "tight", pad_inches = 0.1);
PyPlot.clf();
#--------------------------------------------------
F_tf = abs2.(sum(conj.(opr.wfc_t[100,:]).*opr.wfc_wf).*prt.dx);

writedlm("data/F_tf_"*string(opr.tf)*".txt",F_tf);

open("Fid_sta_tf0.1.txt", "a") do io
       writedlm(io, F_tf)
       end
#-----------------------------end------------------------------------------
