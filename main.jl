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

Read the input from the input file, creates all the global numbering variables used in later stages

Inputs:
	inp_file- Name of the input YAML file with .yaml extension
    debug   - Integer 0,1 to switch on or off debugging

Outputs:
	sref      - reference surface area
	bref      - reference span 
	cref      - reference chord length
	aseq      - sequence of alpha to investigate of form [starting,ending,number of steps] of variable types [float,float,integer]
	vinf      - value of vinf
	ρ         - density
	alpha     - angle of attack
	nwing     - number of wings
	ntot_patch- total number of patches
	npatch_perwing
	          - patches per wing (used for global patches numbering)
	ntot_lat  - total number of lattices (used for global lattice numbering)
	nspan_perpatch
	          - number of spanwise strips per patch (used to loop over spanwise strips)
	nlat_perpatch
			  - number of total lattices per patch
	xyz_qc_patch
			  - quarter chord locations for each patch of form [(x,y,z),(start,end),global patch number]
	chord_patch
		      - chord lengths for the pacthes of form [(start,end),global patch number]
	twist_patch
			  - twists for the pacthces of the form [(start,end),gloabl patch number]		  
	 α_zl_patch
	 		  - zero lift alpha values of the form [(start,end),gloabl patch number]

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
                        nlat_perpatch[counter]    = inp_data[str]["patch"*string(j)]["nlat"]*nspan_perpatch[counter]
						
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
	
	Inputs:
		nwing      - number of wings
		ntot_patch - total number of patches  
		npatch_perwing
				   - number of patches per wing, form : [Integer value]
		ntot_lat   - total numer of lattices	
		nspan_perpatch
				   - number of spanwise strips per patch
		nlat_perpatch
				   - number of lattices per patch (uses global patch numbering)	
		xyz_qc_patch
				   - quarter chord locations for each patch of form [(x,y,z),(start,end),global patch number]
		 chord_patch
				   - chord lengths for the pacthes of form [(start,end),global patch number]
		 twist_patch
				   - twists for the pacthces of the form [(start,end),gloabl patch number]		  
		α_zl_patch
	 			   - zero lift alpha values of the form [(start,end),gloabl patch numb
		debug      - debug flag to enable some debugging options
	
	Output:
		nspn       - number of span wise locations to plt local distributions like Cl
						: If more than one wing, each wings span wise locations are stored in that order
						  Example - Wing + Tail 
						  nspn = span locations on wing + span locations on tail
		spn_map    - mapping for all the lattices to span locations 
						: form - Array{Integers},size "ntot_lat" pointing to position in "spn_loc" array
		spn_loc    - locations along span (y-coordinates)
						: form - Array{Float,1}, size "nspn" 
		sbar       - locations of starting point of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		ebar       - locations of ending point of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		mbar       - locations of middle point of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		nbar       - perpendicular unit vectors of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		cbar       - locations of control points or quarter chord points of horse-shoe lattice for all lattices
						: form - Array{Float,2}, size {3,ntot_lat}
		dbar       - locations of middle points of horse-shoe lattice for all lattices at Treffts Plane
						: form - Array{Float,2}, size {3,ntot_lat}
		tbar       - unit vector for pointing the tail of horse shoe vortices
						: form - Array{Float,2}, size {3,ntot_lat} 
		chord      - chord length values at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		twist      - twist values at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		α_zl       - alpha zero lift values at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		θ          - dihedral angle at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		Λ          - sweep angle at all lattices
						: form - Array{Float,1}, size {ntot_lat}
		ds         - length of filament for all lattices
						: form - Array{Float,1}, size {ntot_lat}

"""
function geom_calc(nwing,ntot_patch,npatch_perwing,ntot_lat,nspan_perpatch,nlat_perpatch,xyz_qc_patch,chord_patch,twist_patch,α_zl_patch,debug=0)

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

		
		spn_map = zeros(Int64,ntot_lat) #Spanwise mapping

		
		nlatspn_perwing = zeros(Int64,nwing) #Span per wing

		# Loop over to find total span wise locations (Could be eliminated/improved)
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

		
		spn_loc = zeros(nspn) #Spanwise locations
		
		#Initialise counters
        ilat  = 0
        ipatch= 0
		nspn  = 0

		
        for i in  tqdm(1:nwing,leave=false)
                for j in 1:npatch_perwing[i]
                        ipatch = ipatch + 1 #for global patch number

						# Interpolate chord, twist and α_zl values
                        itp_chord = LinearInterpolation(xyz_qc_patch[2,:,ipatch],chord_patch[:,ipatch])
                        itp_twist = LinearInterpolation(xyz_qc_patch[2,:,ipatch],twist_patch[:,ipatch])
                        itp_α_zl  = LinearInterpolation(xyz_qc_patch[2,:,ipatch],α_zl_patch[:,ipatch])

						# Get the starting and ending locations along the chord for each patch
						nlat_perspan = Int64(nlat_perpatch[ipatch]/nspan_perpatch[ipatch])
						a = copy(xyz_qc_patch[:,1,ipatch])
                        b = copy(xyz_qc_patch[:,2,ipatch])
						
						#
						# Computing starting root and tip locations for each span
						#
						a_span = zeros(3,nspan_perpatch[ipatch])
						b_span = zeros(3,nspan_perpatch[ipatch])
						
						#pvec is perpendicular to (b-a)
						pvec   = (b[1:2] - a[1:2])/norm(b[1:2]-a[1:2])
						pvec   = [pvec[2],-pvec[1]]
						
						#loop over strips along the chord wise
						ns     = nspan_perpatch[ipatch]
						for ispan in 1:ns
							a_span[:,ispan] = a + 0.25*chord_patch[1,ipatch]*(4*ispan-3-ns)/ns * [pvec[1],pvec[2],0.0]
							b_span[:,ispan] = b + 0.25*chord_patch[2,ipatch]*(4*ispan-3-ns)/ns * [pvec[1],pvec[2],0.0]
						end
						
						#Sanity check to see the pvec vector is in the right direction
						@info pvec
						if (dot(pvec,[1.0,0.0])<=0.0)
							println("Error in geom_calc span wise variables")
							println("Perp Vector Calc wrong")
							return 
						end

						#Loop over the chord wise strips and lattices in each strip					
						for ispan in 1:nspan_perpatch[ipatch]
						for k in 1:nlat_perspan
								
                                ilat = ilat + 1
								
                                a = copy(a_span[:,ispan])
                                b = copy(b_span[:,ispan])
                                
								#Find start,end,mid,normal vec, dihedral, sweep, α_zl, ...
								sbar[:,ilat] =  a + (b - a)/nlat_perspan .*(k-1)
                                ebar[:,ilat] =  a + (b - a)/nlat_perspan .*(k)
                                mbar[:,ilat] = 0.5*(sbar[:,ilat]+ebar[:,ilat])
                                nbar[:,ilat] = cross((ebar[:,ilat] - sbar[:,ilat]),[-1.0,0.0,0.0])/norm(cross((ebar[:,ilat] - sbar[:,ilat]),[-1.0,0.0,0.0]))

                                θ[ilat]      = acos(nbar[3,ilat])
								Λ[ilat]      = acos(dot(b[1:2]-a[1:2],[0.0,1.0])/norm(b[1:2]-a[1:2]))
                                chord[ilat]  = itp_chord(mbar[2,ilat])
                                twist[ilat]  = itp_twist(mbar[2,ilat])
                                α_zl[ilat]   = itp_α_zl(mbar[2,ilat])
								ds[ilat]	 = norm(ebar[2:3,ilat]-sbar[2:3,ilat])#*cos(Λ[ilat])

                                cbar[:,ilat] = mbar[:,ilat] + [chord[ilat]/(2*nspan_perpatch[ipatch]),0.0,0.0]
                                dbar[:,ilat] = [10^8,mbar[2,ilat],mbar[3,ilat]]
								tbar[:,ilat] = [1.0,0.0,0.0]   #Trailing edge vector
								
								#To get the spn_loc locations only when using the first strip
								if (ispan==1)
									nspn = nspn + 1
									spn_loc[nspn] = mbar[2,ilat]
								end

								#Compute ispn_bef to locate where in spn_map this patch is starting
								ispn_bef = 0
								for ip in 1:ipatch-1
									ispn_bef = ispn_bef + Int64(nlat_perpatch[ip]/nspan_perpatch[ip])
								end
								spn_map[ilat] = ispn_bef + k
								
						end
						end
					#finding numer of lattices in spn wise per each wing
					nlatspn_perwing[i] = nlatspn_perwing[i] + nlat_perspan  
				end	
        end

		#Debug stuff
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

    Function to influence for semi-infinte vortex element given starting,ending and tail direction

    Input :
		rbar - vector for pointing the location of evaluation
				: form - Array{Float,1}, size {3} 
		sbar - vector for pointing the start of filament
				: form - Array{Float,1}, size {3} 
		tbar - vector for pointing the trailing of vortex 
				: form - Array{Float,1}, size {3}

    Output:
		q    - vector of influence at the given location
				: form - Array{Float,1}, size {3}

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

    Function to find the three velocity componenets for downwash given the line segment and the point

	Input :
		rbar - vector for pointing the location of evaluation
				: form - Array{Float,1}, size {3} 
		sbar - vector for pointing the start of filament
				: form - Array{Float,1}, size {3} 
		ebar - vector for pointing the end of filament 
				: form - Array{Float,1}, size {3}

    Output:
		q    - vector of influence at the given location
				: form - Array{Float,1}, size {3}	
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
		Γ = vinf*(AIC\rhs)

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


Cl_spn,nspn,spn_map,spn_loc,θ,rhs,AIC,AICₜ,AICₘ,Λ,sbar,ebar,ds,Γ,chord,cbar,mbar,Cl,CL,CDind,CDind_ff,SpnLd = main("input.yaml");


#Plot the local Cl
plt = plot(spn_loc[1:200],Cl_spn[1:200],label="Wing") #wing
plot!(spn_loc[201:260],Cl_spn[201:260],label="Tail")  #tail
plot!(title="Cl vs span",xlabel="Span wise location",ylabel="Cl")
savefig("wingtail_test.png")
display(plt)

#plt = scatter(spn_loc[1:200],Cl_spn[1:200],label="Wing") #wing
#scatter!(spn_loc[201:260],Cl_spn[201:260],label="Tail")  #tail
#scatter!(title="Cl vs span",xlabel="Span wise location",ylabel="Cl")
#savefig("wingtail_test.png")
#display(plt)