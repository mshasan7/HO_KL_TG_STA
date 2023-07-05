function omega_square_STA(prt::Param,opr::Operators,t)
    b(t) = (1.0 .+ 10.0.*(prt.gamma-1).*(t./opr.tf).^3 .- 15.0.*(prt.gamma-1).*(t./opr.tf).^4 .+ 6.0.*(prt.gamma-1).*(t./opr.tf).^5)
    b2(t) = 60.0.*(prt.gamma-1).*t.*(1.0/opr.tf).^3 .- 180.0.*(prt.gamma-1).*t.^2 .*(1.0/opr.tf).^4 .+ 120.0.*(prt.gamma-1).*t.^3 .*(1.0/opr.tf).^5
    w_2(t) = (prt.omega_0.^2 /b(t).^4) .- b2(t)./b(t);
    return w_2(t);
end

function potential_RT_STA(prt::Param,opr::Operators,t) # potential initialisation
    V = 0.5.*m.*omega_square_STA(prt,opr,t).*prt.x.^2;
    return V;
end

function f_transform_RT_STA(prt::Param,opr::Operators,itr)
        forward_fft = plan_fft!(opr.wfc);
        inverse_fft = plan_ifft!(opr.wfc);
        #fft to momentum space
        opr.wfc .= forward_fft*opr.wfc;
        #multiply with the kinetic exponent with half time step
        opr.wfc .= opr.K .*opr.wfc;
        #back to real space
        opr.wfc .= inverse_fft*opr.wfc;
        #potential exponet for the SOC gas
        opr.R .= exp.(-opr.im_value.*(prt.g.*(abs2.(opr.wfc)).^2 .+ potential_RT_STA(prt,opr,(itr-1).*opr.dt)).*opr.dt)
        # Full-step in real space
        opr.wfc .= opr.R .*opr.wfc;
        #fft to momentum space
        opr.wfc .= forward_fft*opr.wfc;
        #multiply with the kinetic exponent with half time step
        opr.wfc .= opr.K .*opr.wfc;
        #back to real space
        opr.wfc .= inverse_fft*opr.wfc;
end

pot_funtion_w0(x) = 0.5 .*m .*prt.omega_0 .^2*(x).^2; # Harmonic trap
pot_funtion_wf(x) = 0.5 .*m .*prt.omega_f.^2*(x).^2; # Harmonic trap with w_f
function potential_w0(prt::Param) # potential initialisation
    V = [pot_funtion_w0(i) for i in prt.x];
    return V;
end
function potential_wf(prt::Param) # potential initialisation
    V = [pot_funtion_wf(i) for i in prt.x];
    return V;
end
function initial_psi(prt::Param,psi_initial::Array{ComplexF64,1})
    @. psi_initial = sqrt(prt.N_bec).*exp(-(prt.x^2));
    norm = sqrt(sum(abs2.(psi_initial)) * prt.dx).*(1/sqrt(prt.N_bec)); #Normalisation
    psi_initial ./= norm;
    return psi_initial;
end
function bec_wavefunction(prt::Param)
    psi = zeros(ComplexF64, prt.nx) ;
    psi = initial_psi(prt, psi);
    return psi;
end
function normalise!(prt::Param,opr::Operators)
    norm = sqrt(sum(abs2.(opr.wfc)).* prt.dx).*(1/sqrt(prt.N_bec));
    opr.wfc ./= norm;
    norm_wf = sqrt(sum(abs2.(opr.wfc_wf)).* prt.dx).*(1/sqrt(prt.N_bec));
    opr.wfc_wf ./= norm_wf;
    return(opr.wfc,opr.wfc_wf)
end


function calculate_energy(prt::Param,opr::Operators,itr)
    energy_final = 0.
    if(imag(opr.im_value) == 0) #energy for imaginary time evolution
        wfc_r = opr.wfc;
        wfc_k = fft(wfc_r);
        wfc_c = conj(wfc_r);
        #kinetic energies
        energy_kin = wfc_c.*ifft((0.5.*(prt.kx.^2)).*wfc_k);
        #potential energies
        energy_pot1 = (1/3).*wfc_c.*(prt.g.*(abs2.(wfc_r)).^2).*wfc_r; #1/3 from the derivation of dE/dphi*;
        energy_pot2 = wfc_c.*(opr.V).*wfc_r;
        energy_pot = energy_pot1 .+ energy_pot2;
        #total energy
        energy = energy_kin .+ energy_pot;
        #energy for V_wf
        wfc_r_wf = opr.wfc_wf;
        wfc_k_wf = fft(wfc_r_wf);
        wfc_c_wf = conj(wfc_r_wf);
        #kinetic energies
        energy_kin_wf = wfc_c_wf.*ifft((0.5.*(prt.kx.^2)).*wfc_k_wf);
        #potential energies
        energy_pot1_wf = (1/3).*wfc_c_wf.*(prt.g.*(abs2.(wfc_r_wf)).^2).*wfc_r_wf; #1/3 from the derivation of dE/dphi*;
        energy_pot2_wf = wfc_c_wf.*(opr.V_wf).*wfc_r_wf;
        energy_pot_wf = energy_pot1_wf .+ energy_pot2_wf;
        #total energy
        energy_wf = energy_kin_wf .+ energy_pot_wf;

        energy_final = real.(sum(energy .+ energy_wf).*prt.dx);
    elseif(imag(opr.im_value) !=0) #calculating energy for real time evolution
        wfc_r = opr.wfc;
        wfc_k = fft(wfc_r);
        wfc_c = conj(wfc_r);
        #kinetic energies
        energy_kin = wfc_c.*ifft((0.5.*(prt.kx.^2)).*wfc_k);
        #potential energies
        energy_pot1 = (1/3).*wfc_c.*(prt.g.*(abs2.(wfc_r)).^2).*wfc_r; #1/3 from the derivation of dE/dphi*;
        energy_pot2 = wfc_c.*(potential_RT_STA(prt,opr,(itr-1).*opr.dt)).*wfc_r;#???????
        energy_pot = energy_pot1 .+ energy_pot2;
        #total energy
        energy = energy_kin .+ energy_pot;
        energy_final = real.(sum(energy).*prt.dx);
    end
    return energy_final;
end

function f_transform_GS(prt::Param,opr::Operators)
        forward_fft = plan_fft!(opr.wfc);
        inverse_fft = plan_ifft!(opr.wfc);
        forward_fft_wf = plan_fft!(opr.wfc_wf);
        inverse_fft_wf = plan_ifft!(opr.wfc_wf);
        #fft to momentum space
        opr.wfc .= forward_fft*opr.wfc;
        opr.wfc_wf .= forward_fft_wf*opr.wfc_wf;
        #multiply with the kinetic exponent with half time step
        opr.wfc .= opr.K .*opr.wfc;
        opr.wfc_wf .= opr.K_wf .*opr.wfc_wf;
        #back to real space
        opr.wfc .= inverse_fft*opr.wfc;
        opr.wfc_wf .= inverse_fft_wf*opr.wfc_wf;
        #potential exponet for the SOC gas
        opr.R .= exp.(-opr.im_value.*(prt.g.*(abs2.(opr.wfc)).^2 .+ opr.V).*opr.dt);
        opr.R_wf .= exp.(-opr.im_value.*(prt.g.*(abs2.(opr.wfc_wf)).^2 .+ opr.V_wf).*opr.dt);
        # Full-step in real space
        opr.wfc .= opr.R .*opr.wfc;
        opr.wfc_wf .= opr.R_wf .*opr.wfc_wf;
        #fft to momentum space
        opr.wfc .= forward_fft*opr.wfc;
        opr.wfc_wf .= forward_fft_wf*opr.wfc_wf;
        #multiply with the kinetic exponent with half time step
        opr.wfc .= opr.K .*opr.wfc;
        opr.wfc_wf .= opr.K_wf .*opr.wfc_wf;
        #back to real space
        opr.wfc .= inverse_fft*opr.wfc;
        opr.wfc_wf .= inverse_fft_wf*opr.wfc_wf;
end

function split_method!(prt::Param, opr::Operators)
    iterations_tot = floor(Int64, opr.tf/opr.dt); #total time for iterations in real time
    jj=1;# counter for saving the wavefunction as a function of time
    if (imag(opr.im_value) == 0)
        iterations_tot = Int64(5e6); # total iteration for imaginary time
        E_temp = 0;
        println("calculating GS...");
    else
        println("Running time evolution...");
    end
    for itr = 1:iterations_tot
        if (imag(opr.im_value) == 0)
            f_transform_GS(prt,opr) #fourier transform routine
            normalise!(prt,opr) # renormalizing for imaginary time
            if (itr % 300 == 0)
                current_E = calculate_energy(prt,opr,itr);  #calculate_energy consumes a lot of memory
                dE = abs(1.0 - abs(E_temp/current_E));
                @printf("energy = %e\t dE = %e\n", E_temp,dE);
                E_temp = current_E;
                if(prt.tolerance > dE ) # Tolerance to stop the program in imaginary time
                    println("Final dE = ", dE);
                    println("Num iter = ", itr);
                    println("Ground state reached");
                    break;
                elseif(itr == iterations_tot)
                    println("Ground state could not be reached");
                    break;
                end
            end
        elseif(imag(opr.im_value) != 0)
                f_transform_RT_STA(prt,opr,itr)
		E_current = calculate_energy(prt,opr,itr) # energy check for the dynamics
                if(mod(itr,(iterations_tot/prt.N_pic))==0)
                    opr.wfc_t[jj,:] .= opr.wfc; #saving the wfc as a function of tiem
                    jj += 1;
                end
            	if (itr % 2000 == 0)
                	print("Energy =", E_current,"\t");
            	end
        end
        if(itr==iterations_tot)
            println("Total iter = ", itr);
        end
    end
end
function gs_plots!(prt::Param,opr::Operators)
figure(figsize=(10, 3))
subplot(1,2,1)
PyPlot.plot(prt.x,abs2.(opr.wfc),label="ground state|w0")
PyPlot.plot(prt.x,abs2.(opr.wfc_in),"-.",label="initial")
ylabel("density")
xlabel("position")
legend()
subplot(1,2,2)
PyPlot.plot(prt.x,abs2.(opr.wfc_wf),label="ground state|wf")
PyPlot.plot(prt.x,abs2.(opr.wfc_in),"-.",label="initial")
ylabel("density")
xlabel("position")
legend()
end

function rt_plots!(prt::Param,opr::Operators)
PyPlot.plot(prt.x,abs2.(opr.wfc_t[100,:]),label=L"wfc|_{t=tf}")
PyPlot.plot(prt.x,abs2.(opr.wfc_t[1,:]),"--",label=L"wfc|_{t=0}")
legend()
end
