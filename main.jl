#Libraries used
using YAML
using Parameters
using Interpolations
using LinearAlgebra
using ProgressBars
using Plots
using NPZ

"""
    read_inp(debug=0)

Read the input from the input file

Inputs:

        debug - Integer 0,1 to switch on or off debugging

Outputs:
"""
function read_inp(inp_file,debug=0)
        inp_data = YAML.load_file(inp_file)
         @info "Reading Input File"
        name = inp_data["name"]
        sref = inp_data["sref"]
        bref = inp_data["bref"]
        cref = inp_data["cref"]
        nwing= inp_data["nwing"]

        #create arrays required to store
        npatch_perwing = zeros(Int64,nwing)
         for i in 1:nwing
                str = "wing"*string(i)
                npatch_perwing[i] = inp_data[str]["npatch"]
        end

        ntot_patch = sum(npatch_perwing)
        #get some arrays to store root and tip info
        xyz_qc_patch = zeros(3,2,ntot_patch)
        chord_patch  = zeros(2,ntot_patch)
        twist_patch  = zeros(2,ntot_patch)
        α_zl_patch   = zeros(2,ntot_patch)

        nlat_perpatch = zeros(Int64,ntot_patch)
		nspan_perpatch = zeros(Int64,ntot_patch)
        ntot_lat   = 0
        counter    = 0
        for i in 1:nwing
               str = "wing"*string(i)
                for j in 1:npatch_perwing[i]
                        counter = counter + 1
						nspan_perpatch[counter]   = inp_data[str]["patch"*string(j)]["nspan"]
                        nlat_perpatch[counter]    = inp_data[str]["patch"*string(j)]["nlat"]*nspan_perpatch[j]
						
                        # at root
                        xyz_qc_patch[:,1,counter] = inp_data[str]["patch"*string(j)]["root"]["xyz_qc"]
                        chord_patch[1,counter]    = inp_data[str]["patch"*string(j)]["root"]["chord"]
                        twist_patch[1,counter]    = inp_data[str]["patch"*string(j)]["root"]["twist"]*pi/180
                        α_zl_patch[1,counter]     = inp_data[str]["patch"*string(j)]["root"]["alpha_zl"]*pi/180

                        # at tip
                        xyz_qc_patch[:,2,counter] = inp_data[str]["patch"*string(j)]["tip"]["xyz_qc"]
                        chord_patch[2,counter]    = inp_data[str]["patch"*string(j)]["tip"]["chord"]
                        twist_patch[2,counter]    = inp_data[str]["patch"*string(j)]["tip"]["twist"]*pi/180
                        α_zl_patch[2,counter]     = inp_data[str]["patch"*string(j)]["tip"]["alpha_zl"]*pi/180

                        ntot_lat = ntot_lat + nlat_perpatch[counter]
                end
        end

		aseq = inp_data["aseq"]
		aseq = [aseq[1]*pi/180,aseq[2]*pi/180,aseq[3]]
		vinf = inp_data["vinf"]
		ρ    = inp_data["rho"]
		alpha= inp_data["alpha"]*pi/180


        if (debug==1)
                #checks

                #print some stuff if needed
                display("ntot_lat:"*string(ntot_lat))
                display(xyz_qc_patch)

        end

        return (sref,bref,cref,aseq,vinf,ρ,alpha,nwing,ntot_patch,npatch_perwing,ntot_lat,nspan_perpatch,nlat_perpatch,xyz_qc_patch,chord_patch,twist_patch,α_zl_patch)
end

"""
    geom_calc(nwing,npatch,npatch_perwing,ntot_lat,nlat_perpatch,xyz_qc_patch,chord_patch,twist_patch,α_zl_patch,Λ,debug)

Calculate the geometric propoerties for all the lattice elements
"""
function geom_calc(nwing,npatch,npatch_perwing,ntot_lat,nspan_perpatch,nlat_perpatch,xyz_qc_patch,chord_patch,twist_patch,α_zl_patch,debug=0)

        #create arrays for the lattices
        sbar = zeros(3,ntot_lat)        # starting of bv
        ebar = zeros(3,ntot_lat)        # ending of bv
        mbar = zeros(3,ntot_lat)        # mid of bv
        cbar = zeros(3,ntot_lat)        # ctrl poi of bv
        dbar = zeros(3,ntot_lat)        # trefftz plane mid poi
        nbar = zeros(3,ntot_lat)		# normal for patch/wing/lattice
        tbar = zeros(3,ntot_lat)        # tail direction

        chord = zeros(ntot_lat)         # chord length
        twist = zeros(ntot_lat)			# twist
        α_zl  = zeros(ntot_lat)			# α_zero-lift
        θ     = zeros(ntot_lat)			# dihedral angle
		Λ     = zeros(ntot_lat)         # sweep angle
		ds    = zeros(ntot_lat)         # lattice length

		#Spanwise mapping
		spn_map = zeros(Int64,ntot_lat)

		#Span per wing
		nlatspn_perwing = zeros(Int64,nwing)

		nspn  = 0
		ipatch=0
		for i in 1:nwing
			for j in 1:npatch_perwing[i]
				ipatch = ipatch + 1
				nlat_perspan = Int64(nlat_perpatch[ipatch]/nspan_perpatch[ipatch])
				for k in 1:nlat_perspan
					nspn = nspn + 1			
				end
			end
		end

		#Spanwise locations
		spn_loc = zeros(nspn)

        ilat  = 0
        ipatch= 0
		nspn  = 0
		#Trailing edge
		tbar[1:3,:] .= [1.0,0.0,0.0]
        for i in  tqdm(1:nwing,leave=false)
                for j in 1:npatch_perwing[i]
                        ipatch = ipatch + 1
                        itp_chord = LinearInterpolation(xyz_qc_patch[2,:,ipatch],chord_patch[:,ipatch])
                        itp_twist = LinearInterpolation(xyz_qc_patch[2,:,ipatch],twist_patch[:,ipatch])
                        itp_α     = LinearInterpolation(xyz_qc_patch[2,:,ipatch],α_zl_patch[:,ipatch])

						nlat_perspan = Int64(nlat_perpatch[ipatch]/nspan_perpatch[ipatch])
						a = copy(xyz_qc_patch[:,1,ipatch])
                        b = copy(xyz_qc_patch[:,2,ipatch])
						a_span = zeros(3,nspan_perpatch[ipatch])
						b_span = zeros(3,nspan_perpatch[ipatch])
						pvec   = (b[1:2] - a[1:2])/norm(b[1:2]-a[1:2])
						pvec   = [pvec[2],-pvec[1]]
						ns     = nspan_perpatch[ipatch]
						for ispan in 1:ns
							a_span[:,ispan] = a + 0.25*chord_patch[1,ipatch]*(4*ispan-3-ns)/ns * [pvec[1],pvec[2],0.0]
							b_span[:,ispan] = b + 0.25*chord_patch[2,ipatch]*(4*ispan-3-ns)/ns * [pvec[1],pvec[2],0.0]
						end
						#Sanity check
						#@info pvec
						if (dot(pvec,[1.0,0.0])<=0.0)
							println("Error in geom_calc span wise variables")
							println("Perp Vector Calc wrong")
							return 
						end
						#@info 1,nspan_perpatch[:]
						for ispan in 1:nspan_perpatch[ipatch]
                        for k in 1:nlat_perspan
                                ilat = ilat + 1
                                a = copy(a_span[:,ispan])
                                b = copy(b_span[:,ispan])
                                sbar[:,ilat] =  a + (b - a)/nlat_perspan .*(k-1)
                                ebar[:,ilat] =  a + (b - a)/nlat_perspan .*(k)
                                mbar[:,ilat] = 0.5*(sbar[:,ilat]+ebar[:,ilat])
                                nbar[:,ilat] = cross((ebar[:,ilat] - sbar[:,ilat]),[-1.0,0.0,0.0])/norm(cross((ebar[:,ilat] - sbar[:,ilat]),[-1.0,0.0,0.0]))

                                θ[ilat]      = acos(nbar[3,ilat])
								Λ[ilat]      = acos(dot(b[1:2]-a[1:2],[0.0,1.0])/norm(b[1:2]-a[1:2]))
                                chord[ilat]  = itp_chord(mbar[2,ilat])
                                twist[ilat]  = itp_twist(mbar[2,ilat])
                                α_zl[ilat]   = itp_α(mbar[2,ilat])
								ds[ilat]	 = norm(ebar[2:3,ilat]-sbar[2:3,ilat])#*cos(Λ[ilat])

                                cbar[:,ilat] = mbar[:,ilat] + [chord[ilat]/2,0.0,0.0]
                                dbar[:,ilat] = [10^8,mbar[2,ilat],mbar[3,ilat]]

								if (ispan==1)
									nspn = nspn + 1
									#@info 2,nspn
									spn_loc[nspn] = mbar[2,ilat]

								end

								#@info i,ispan,j,k
								ilat_bef = 0
								for ip in 1:ipatch-1
									ilat_bef = ilat_bef + Int64(nlat_perpatch[ip]/nspan_perpatch[ip])
								end
								spn_map[ilat] = ilat_bef + k
								
						end
						end
					nlatspn_perwing[i] = nlatspn_perwing[i] + nlat_perspan  
				end	
        end
		#@info nlatspn_perwing

        if debug==1
                display("sbar")
                display(sbar)
                display("xyz_qc_patch")
                display(xyz_qc_patch)
        end
        return nspn,spn_map,spn_loc,sbar,ebar,mbar,nbar,cbar,dbar,tbar,chord,twist,α_zl,θ,Λ,ds
end
"""
    calc_vind_sinf(rbar,sbar,tbar)

    Input :

    Output:

    Function to influence for semi-infinte vortex element given starting,ending and tail direction
"""
function calc_q_sinf(rbar,sbar,tbar)

	(x1,y1,z1) = sbar
	(x,y,z)    = rbar

	avec = [x1-x,y1-y,z1-z]


	amag = norm(avec)


	cprd = cross(avec,tbar)
	dprd = dot(avec,tbar)

	num = cprd* (1.0-(dprd/(amag)))
	den = dot(cprd,cprd)

	#for the point coinciding with the element
	if (den<1e-6)
		q = [0.0,0.0,0.0]
	else
		q =  num/den
	end
	return q
end
"""
    calc_vind_finite(rbar,sbar,ebar)

    Input :

    Output:

    Function to find the three velocity componenets for downwash given the line segment and the point
"""
function calc_q_finite(rbar,sbar,ebar)

	(x1,y1,z1) = sbar
	(x2,y2,z2) = ebar
	(x,y,z)    = rbar

	avec = [x1-x,y1-y,z1-z]
	bvec = [x2-x,y2-y,z2-z]

	amag = norm(avec)
	bmag = norm(bvec)

	cprd = cross(avec,bvec)
	dprd = dot(avec,bvec)

	num = cprd* (amag+bmag) * (1.0-(dprd/(amag*bmag)))
	den = dot(cprd,cprd)

	#for the point coinciding with the element
	if (den<1e-6)
		q = [0.0,0.0,0.0]
	else
		q =  num/den
	end
	return q
end
"""
    calc_AICs()

    Calculate the AIC matrices required
"""
function calc_AICs(ntot_lat,sbar,ebar,cbar,mbar,nbar,tbar,debug=0)

	pi_inv = 1.0/pi

        AIC  = zeros(ntot_lat,ntot_lat)
        AICₘ = zeros(ntot_lat,ntot_lat)
        AICₜ = zeros(ntot_lat,ntot_lat)

	for i in 1:ntot_lat
		for j in 1:ntot_lat

			#Calculating AIC
			#bounded vortex contribution
			qbv     = calc_q_finite(cbar[:,i],sbar[:,j],ebar[:,j])
			qbv_comp= dot(qbv,nbar[:,i])
			#left leg contribution
			qlv     = -calc_q_sinf(cbar[:,i],sbar[:,j],tbar[:,j])
			qlv_comp= dot(qlv,nbar[:,i])
			#right leg contribution
			qrv     = calc_q_sinf(cbar[:,i],ebar[:,j],tbar[:,j])
			qrv_comp= dot(qrv,nbar[:,i])

			AIC[i,j] = 0.25*pi_inv*(qbv_comp+qlv_comp+qrv_comp)
			qbv .= 0.0
			qbv_comp = 0.0
			qlv .= 0.0
			qlv_comp = 0.0
			qrv .= 0.0
			qrv_comp = 0.0

			#Calculating AICₘ
			#bounded vortex contribution
			qbv     = calc_q_finite(mbar[:,i],sbar[:,j],ebar[:,j])
			qbv_comp= dot(qbv,-nbar[:,i])
			#left leg contribution
			qlv     = -calc_q_sinf(mbar[:,i],sbar[:,j],tbar[:,j])
			qlv_comp= dot(qlv,-nbar[:,i])
			#right leg contribution
			qrv     = calc_q_sinf(mbar[:,i],ebar[:,j],tbar[:,j])
			qrv_comp= dot(qrv,-nbar[:,i])

			AICₘ[i,j] = 0.25*pi_inv*(qbv_comp+qlv_comp+qrv_comp)

			qbv .= 0.0
			qbv_comp = 0.0
			qlv .= 0.0
			qlv_comp = 0.0
			qrv .= 0.0
			qrv_comp = 0.0

			#Calculating AICₜ
			#coordinate transformation for calculation purposes
			c0 = cbar[:,i]
			s0 = sbar[:,j]
			e0 = ebar[:,j]

			c0[1] = s0[1] = e0[1] = 0.0

			#bounded vortex contribution
			qbv     = 0.0
			#left leg contribution
			qlv     = -calc_q_sinf(c0,s0,tbar[:,j])
			qlv_comp= dot(qlv,nbar[:,i])
			#right leg contribution
			qrv     = calc_q_sinf(c0,e0,tbar[:,j])
			qrv_comp= dot(qrv,nbar[:,i])

			AICₜ[i,j] = 0.25*pi_inv*(qbv_comp+2.0*qlv_comp+2.0*qrv_comp)
		end
	end


        return AIC,AICₘ,AICₜ
end
"""
    calc_rhs()

	Function to calculate the rhs

	Input :

	Output:

"""
function calc_rhs(ntot_lat,α,θ,twist,α_0l)
	rhs = zeros(ntot_lat)
	for ilat in 1:ntot_lat
		α_g	 = asin(sin(α)*cos(θ[ilat]))
		rhs[ilat] = -sin(α_g+twist[ilat]-α_0l[ilat])
	end

	return rhs
end
"""
    main(inp_file)

	Main function to calculate all the values required

	Input : inp_file - Input File

	Output:
"""
function main(inp_file,iseq=0)

    sref,bref,cref,aseq,vinf,ρ,alpha,nwing,
	ntot_patch,npatch_perwing,
	ntot_lat,nspan_perpatch,nlat_perpatch,xyz_qc_patch,
	chord_patch,twist_patch,α_zl_patch = read_inp(inp_file)

	nspn,spn_map,spn_loc,sbar,ebar,mbar,nbar,cbar,dbar,tbar,chord,twist,α_zl,θ,Λ,ds = geom_calc(nwing,ntot_patch,npatch_perwing,ntot_lat,nspan_perpatch,nlat_perpatch,xyz_qc_patch,chord_patch,twist_patch,α_zl_patch)

	AIC,AICₘ,AICₜ = calc_AICs(ntot_lat,sbar,ebar,cbar,mbar,nbar,tbar)

	if iseq == 1
		α = LinRange(aseq[1],aseq[2],Int64(aseq[3]))
		nα = size(α,1)

		CL_seq       = zeros(nα)
		CDind_seq    = zeros(nα)
		CDind_ff_seq = zeros(nα)

 		for i in 1:nα
			rhs = calc_rhs(ntot_lat,α[i],θ,twist,α_zl)
			Γ = vinf*AIC\rhs

			# Calculate lift and lift coefficient
			dL = 0.5*ρ*vinf*Γ.*ds # only for static
			Cl = 2.0*dL./(ρ*(vinf^2)*chord.*ds)
			CL = 2.0/sref * sum(Γ.*cos.(θ).*ds/(vinf))

			# Calculate induced drag and induced drag coefficient
			w = AICₘ*Γ
			Dind = ρ*sum(w.*Γ.*ds)
			CDind= 2.0/(vinf^2*sref)*sum(w.*Γ.*ds)

			# Calculate induced drag and coefficient at far field
			vₜ = -AICₜ * Γ
			Dind_ff = 0.5*ρ*sum(Γ.*vₜ.*ds)
			CDind_ff= 1.0/(vinf^2*sref)*sum(vₜ.*Γ.*ds)

			#Calculate Span load
			cavg  = sum(chord.*ds)/(sum(ds))
			SpnLd = chord.*Cl./(cavg*CL)

			CL_seq[i] = CL
			CDind_seq[i] = CDind
			CDind_ff_seq[i] = CDind_ff

		end

		#calculate slope
		slp_CL     = (CL_seq[nα] - CL_seq[1])/(α[nα] - α[1])
		slp_CDind   = (CDind_seq[nα] - CDind_seq[1])/(α[nα] - α[1])
		slp_CDind_ff= (CDind_ff_seq[nα] - CDind_ff_seq[1])/(α[nα] - α[1])

		return α,CL_seq,CDind_seq,CDind_ff_seq,slp_CL,slp_CDind,slp_CDind_ff

	else
		rhs = calc_rhs(ntot_lat,alpha,θ,twist,α_zl)
		Γ = vinf*AIC\rhs

		# Calculate lift and lift coefficient
		dL = 0.5*ρ*vinf*Γ.*ds # only for static
		Cl = 2.0*dL./(ρ*(vinf^2)*chord.*ds)
		CL = 2.0/sref * sum(Γ.*cos.(θ).*ds/(vinf))

		# Calculate induced drag and induced drag coefficient
		w = AICₘ*Γ
		Dind = ρ*sum(w.*Γ.*ds)
		CDind= 2.0/(vinf^2*sref)*sum(w.*Γ.*ds)

		# Calculate induced drag and coefficient at far field
		vₜ = -AICₜ * Γ
		Dind_ff = 0.5*ρ*sum(Γ.*vₜ.*ds)
		CDind_ff= 1.0/(vinf^2*sref)*sum(vₜ.*Γ.*ds)

		#Calculate Span load
		cavg  = sum(chord.*ds)/(sum(ds))
		SpnLd = chord.*Cl./(cavg*CL)

		#Calculate values over span
		Cl_spn = zeros(nspn)
		for ilat in 1:ntot_lat
			ispn = spn_map[ilat]
			#@info ispn
			Cl_spn[ispn] = Cl_spn[ispn] + Cl[ilat] 
		end

		return Cl_spn,nspn,spn_map,spn_loc,θ,rhs,AIC,AICₜ,AICₘ,Λ,sbar,ebar,ds,Γ,chord,cbar,mbar,Cl,CL,CDind,CDind_ff,SpnLd
	end
end

#%%

#cbar,Cl,CL,CDind,CDind_ff,SpnLd = main("input.yaml")
#α,CL_seq,CDind_seq,CDind_ff_seq,slp_CL,slp_CDind,slp_CDind_ff = main("input.yaml",1)
Cl_spn,nspn,spn_map,spn_loc,θ,rhs,AIC,AICₜ,AICₘ,Λ,sbar,ebar,ds,Γ,chord,cbar,mbar,Cl,CL,CDind,CDind_ff,SpnLd = main("input.yaml");

 #%%
 #plot(xlabel="x-coord",ylabel="y-coord",zlabel="Cl-Local Lift")
 #plot!(mbar[1,:],mbar[2,:],Cl,st=:surface,label="Cl")
 #plot!(mbar[1,:],mbar[2,:],Cl,color=:black,camera=(60,30),label="Local Lift Coeff")

 plot(spn_loc[1:100],Cl_spn[1:100])
 plot!(spn_loc[101:154],Cl_spn[101:154])
 savefig("wingtail_test.png")