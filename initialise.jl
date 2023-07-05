struct Param
    nx:: Int64
    dx:: Float64
    omega_0:: Float64
    omega_f:: Float64
    gamma:: Float64
    g:: Float64
    N_bec:: Int64
    tolerance:: Float64
    totspace:: Float64
    x:: Vector{Float64}
    dk:: Float64
    kx:: Array{Float64}
    im_time:: Bool
    re_time:: Bool
    img_output:: Bool
    save_wfc:: Bool
    N_pic:: Int64
    Param(my_dict::Dict{String,Any}) = new(
    my_dict["nx"],
    my_dict["dx"],
    my_dict["omega_0"],
    my_dict["omega_f"],
    my_dict["gamma"],
    my_dict["g"],
    my_dict["N_bec"],
    my_dict["tolerance"],
    my_dict["totspace"],
    collect(-my_dict["totspace"]/2.0+my_dict["dx"]/2:my_dict["dx"]:my_dict["totspace"]/2.0),
    pi/(my_dict["totspace"]/2),
    Vector{Float64}(vcat(0:my_dict["nx"]/2-1, -my_dict["nx"]/2:-1)*pi/(my_dict["totspace"]/2.0)),
    my_dict["imaginary_time"],
    my_dict["real_time"],
    my_dict["images_output"],
    my_dict["save_wfc"],
    my_dict["N_pic"]
    )
    end


mutable struct Operators
    im_value::Complex{Float64}
    V::Array{Float64}
    V_wf::Array{Float64}
    V_RT_STA::Array{Float64}
    R::Array{Complex{Float64}}
    K::Array{Complex{Float64}}
    R_wf::Array{Complex{Float64}}
    K_wf::Array{Complex{Float64}}
    wfc::Array{Complex{Float64}}
    wfc_wf::Array{Complex{Float64}}
    wfc_in::Array{Complex{Float64}}
    wfc_t::Array{Complex{Float64}}
    dt::Float64
    tf::Float64
    Operators(my_dict::Dict{String,Any}) =
        new(0*im,
        zeros(my_dict["nx"]),
        zeros(my_dict["nx"]),
        zeros(my_dict["nx"]),
        zeros(my_dict["nx"]),
        zeros(my_dict["nx"]),
        zeros(my_dict["nx"]),
        zeros(my_dict["nx"]),
        zeros(my_dict["nx"]),
        zeros(my_dict["nx"]),
        zeros(my_dict["nx"]),
        zeros(my_dict["N_pic"],my_dict["nx"]),
        5e-3,
        1.0)
end
