
using ModelingToolkit, DifferentialEquations, Plots


@parameters t
D = Differential(t)


parspop = @parameters r₁=0.038 Kp₁=1.5 r₂=0.042 Kp₂=9.7 mₒ=5 m₁=5 #population/migration dynamics
parsecn = @parameters a = 0.04 b = 0.04 c = 0.01 γ = 0.0016 Ypc0 = 10 δ₁=0.05 δ₂=0.05 α₁=0.5 α₂=0.5 s₁=0.25 s₂=0.21 Cₘ₁=0.7 Cₘ₂=0.7 a₁=2.7 a₂=1.7 #economic dynamics
parsear = @parameters r = 0.1 re₁=0.1 re₂=0.1 eb=0.00004 u=0.0025 α=0.1 σ₁=0.03 σ₂=0.03 G₀=20 G₁=5 Tg=50 η=1 #earth system and climate 
pars = [parspop;parsecn;parsear] #combine into a single vector.

vars = @variables P₁(t)=0.24 P₂(t)=0.24 K₁(t)=0 K₂(t)=0 G(t)=2.8 e₁(t)=0.0004 e₂(t)=0.0004 z(t)=0 P1_nomig(t)=0.24 P2_nomig(t)=0.24 K1_nomig(t)=0 K2_nomig(t)=0 G_nomig(t)=2.8 e1_nomig(t)=0.0004 e2_nomig(t)=0.0004 z_nomig(t)=0 D12(t)=0 D21(t)=0 P₁_gp(t)=0.24 P₂_gp(t)=0.24 K₁_gp(t)=0 K₂_gp(t)=0 G_gp(t)=2.8 e₁_gp(t)=0.0004 e₂_gp(t)=0.0004 z_gp(t)=0 P₁_gpc(t)=0.24 P₂_gpc(t)=0.24 K₁_gpc(t)=0 K₂_gpc(t)=0 G_gpc(t)=2.8 e₁_gpc(t)=0.0004 e₂_gpc(t)=0.0004 z_gpc(t)=0 P₁_gpcr(t)=0.24 P₂_gpcr(t)=0.24 P_kp(t)=0.24 K₁_gpcr(t)=0 K₂_gpcr(t)=0 G_gpcr(t)=2.8 e₁_gpcr(t)=0.0004 e₂_gpcr(t)=0.0004 z_gpcr(t)=0

θ(x) = (max(0,x))^5/(0.0001^5 + (max(0,x))^5)
@register θ(x)

ϕ(x) = x*θ(x)
@register ϕ(x)


Yᵢ(a,α,K,P) = a*((ϕ(K))^α)*((ϕ(P))^(1-α))
@register Yᵢ(a,α,K,P)

Yᵢ₁ = Yᵢ(a₁,α₁,K₁,P₁)
Yᵢ₂ = Yᵢ(a₂,α₂,K₂,P₂)
Yᵢ₁_gp = Yᵢ(a₁,α₁,K₁_gp,P₁_gp)
Yᵢ₂_gp = Yᵢ(a₂,α₂,K₂_gp,P₂_gp)
Yᵢ₁_gpc = Yᵢ(a₁,α₁,K₁_gpc,P₁_gpc)
Yᵢ₂_gpc = Yᵢ(a₂,α₂,K₂_gpc,P₂_gpc)
Yᵢ₁_gpcr = Yᵢ(a₁,α₁,K₁_gpcr,P₁_gpcr)
Yᵢ₂_gpcr = Yᵢ(a₂,α₂,K₂_gpcr,P₂_gpcr)
Yᵢ₁_nomig = Yᵢ(a₁,α₁,K1_nomig,P1_nomig)
Yᵢ₂_nomig = Yᵢ(a₂,α₂,K2_nomig,P2_nomig)

Y(P,Y) = P + ϕ(Y-P)
@register Y(P,K)

Y₁ = Y(P₁,Yᵢ₁)
Y₂ = Y(P₂,Yᵢ₂)
Y₁_gp = Y(P₁_gp,Yᵢ₁_gp)
Y₂_gp = Y(P₂_gp,Yᵢ₂_gp)
Y₁_gpc = Y(P₁_gpc,Yᵢ₁_gpc)
Y₂_gpc = Y(P₂_gpc,Yᵢ₂_gpc)
Y₁_gpcr = Y(P₁_gpcr,Yᵢ₁_gpcr)
Y₂_gpcr = Y(P₂_gpcr,Yᵢ₂_gpcr)
Y₁_nomig = Y(P1_nomig,Yᵢ₁_nomig)
Y₂_nomig = Y(P2_nomig,Yᵢ₂_nomig)

er₁ = θ(P₁-Yᵢ₁)*eb + θ(Yᵢ₁-P₁)*e₁
er₂ = θ(P₂-Yᵢ₂)*eb + θ(Yᵢ₂-P₂)*e₂
er₁_gp = θ(P₁_gp-Yᵢ₁_gp)*eb + θ(Yᵢ₁_gp-P₁_gp)*e₁_gp
er₂_gp = θ(P₂_gp-Yᵢ₂_gp)*eb + θ(Yᵢ₂_gp-P₂_gp)*e₂_gp
er₁_gpc = θ(P₁_gpc-Yᵢ₁_gpc)*eb + θ(Yᵢ₁_gpc-P₁_gpc)*e₁_gpc
er₂_gpc = θ(P₂_gpc-Yᵢ₂_gpc)*eb + θ(Yᵢ₂_gpc-P₂_gpc)*e₂_gpc
er₁_gpcr = θ(P₁_gpcr-Yᵢ₁_gpcr)*eb + θ(Yᵢ₁_gpcr-P₁_gpcr)*e₁_gpcr
er₂_gpcr = θ(P₂_gpcr-Yᵢ₂_gpcr)*eb + θ(Yᵢ₂_gpcr-P₂_gpcr)*e₂_gpcr
er₁_nomig = θ(P1_nomig-Yᵢ₁_nomig)*eb + θ(Yᵢ₁_nomig-P1_nomig)*e1_nomig
er₂_nomig = θ(P2_nomig-Yᵢ₂_nomig)*eb + θ(Yᵢ₂_nomig-P2_nomig)*e2_nomig


#Diaspora(P,M) for remittances 
dia(P,P_nm)=P-P_nm
@register dia(P,P_nm)
D_21_g = 0
D_21_gpc = 0#dia(P₁_gpc,P1_nomig)
D_21_gpcr = 0#dia(P₁_gpcr,P1_nomig)

Yd₁ = ϕ(Y₁ - P₁*Cₘ₁ - r*(D_21_g*(Y₁/P₁)))
Yd₂ = ϕ(Y₂ - P₂*Cₘ₂ + r*(D_21_g*(Y₁/P₁)))
Yd₁_gp = ϕ(Y₁_gp - P₁_gp*Cₘ₁)
Yd₂_gp = ϕ(Y₂_gp - P₂_gp*Cₘ₂)
Yd₁_gpc = ϕ(Y₁_gpc - P₁_gpc*Cₘ₁ - r*(D_21_gpc*(Y₁_gpc/P₁_gpc)))
Yd₂_gpc = ϕ(Y₂_gpc - P₂_gpc*Cₘ₂ + r*(D_21_gpc*(Y₁_gpc/P₁_gpc)))
Yd₁_gpcr = ϕ(Y₁_gpcr - P₁_gpcr*Cₘ₁)
Yd₂_gpcr = ϕ(Y₂_gpcr - P₂_gpcr*Cₘ₂)
Yd₁_nomig = ϕ(Y₁_nomig - P1_nomig*Cₘ₁)
Yd₂_nomig = ϕ(Y₂_nomig - P2_nomig*Cₘ₂)



################ Migration ##################

#weighting factors for
## gdp per capita difference + a
## population size - b
## obstacles - migallow 
## integration +/- i
## diaspora size d
## return migration? 
## remittances r

##population size --> holds back migration  ## argue from literature 
# share of population from region /flobal 
pop_i(Pi,Pj) = 1-(Pi/(Pj+Pi))#Pi./(Pi+Pj)
@register pop_i(Pi,Pj)

pop2_gp = pop_i(P₂_gp,P₁_gp)
pop2_gpc = pop_i(P₂_gpc,P₁_gpc)
pop2_gpcr = pop_i(P₂_gpcr,P₁_gpcr)

#GDP difference function (share of per capita at home vs per capita destination)
#gdp_ij(Yj,Pj,Yi,Pi) = (Yj./(Pj.*((Yj./Pj)+(Yi./Pi))))


## from more people to poor to move
gdp2_ij(Yj,Pj,Yi,Pi) = 1- ((Yi./Pi)/(Yj./Pj))#1/(Yi./(Pi.*((Yj./Pj)+(Yi./Pi))))#. 1/*(1/(1+exp.(γ*((Yi./Pi)-Ypc0))))#
@register gdp2_ij(Yj,Pj,Yi,Pi)

#pij(Yj,Pj,Yi,Pi) = 1 .-(Yi./(Pi.*((Yi./Pi)+(Yj./Pj)))) 
#pij(Yj,Pj,Yi,Pi) = 1 -(Yj/(Pj*((Yj/Pj)+(Yi/Pi)))) ## falschrum
#@register pij(Yj,Pj,Yi,Pi)


## Migrationr rate 
Mij(P,m) = P.*m
@register Mij(P,m)

##depr Carrying capacity for 
kp(Kp,G) = Kp - G*θ(G-4)#(ϕ(G-4))
@register kp(Kp,G)
kp_2 = kp(Kp₂,G_gpc)
kp_2r = kp(Kp₂,G_gpcr)

## migratipon rates 
m21 = gdp2_ij(Y₁,P₁,Y₂,P₂).*a #- pop2.*b + 0.01*θ(P₂-kp_2)
m21_gp = gdp2_ij(Y₁_gp,P₁_gp,Y₂_gp,P₂_gp).*a - pop2_gp.*b
m21_gpc = gdp2_ij(Y₁_gpc,P₁_gpc,Y₂_gpc,P₂_gpc).*a - pop2_gpc.*b #+ c*θ(P_kp-kp_2)
m21_gpcr = gdp2_ij(Y₁_gpcr,P₁_gpcr,Y₂_gpcr,P₂_gpcr).*a - pop2_gpcr.*b #+ 0.01*θ(P_kp-kp_2r)

M21= Mij(P₂,m21)
M21_gp= Mij(P₂_gp,m21_gp)
M21_gpc= Mij(P₂_gpc,m21_gpc)
M21_gpcr= Mij(P₂_gpcr,m21_gpcr)
#define return migration different as push migration? 

##Remittances 
# a share of the Kapital generated in region 1 will be sent 
# to region 2 --> migrants calcualated eacht ime step M12 and M21
# stock of migrants is sum of M12 over time steps 
# M12 is flow, what is stock?
# Assumption: only migrant flows send remittances ?! (assupmtion wrong) 
#r*(M12*(Y₁/P₁)) certain share of per capita income send back home
# idea calculate without migration and substract (stupid)



#migallow, allows migrants to enetr or not 
#--> people disappear or die on the way or get stranded in between 
# --> system looses population
migallow = 1

#differential equations

eqs = [
    ##P1 - only gdp difference
    D(P₁) ~ r₁*P₁*(1 - P₁/Kp₁) + M21*migallow, #human population dynamics 
    D(P₂) ~ r₂*P₂*(1 - P₂/Kp₂) - M21, # standard logistic.
    D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - σ₁*exp(G-G₁)*K₁,# - r*(D_21*(Y₁/P₁)), #change in capital stock = 
    D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - σ₂*exp(G-G₁)*K₂,# + r*(D_21*(Y₁/P₁)), # savings-depreciation-damages
    D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8) + α*θ(G-G₀), # global externality.
    D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
    D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
    D(e₂) ~ -z*re₂*e₂, # decarbonization in region 2
    ### P1_gp gdp and pop 
    D(P₁_gp) ~ r₁*P₁_gp*(1 - P₁_gp/Kp₁) + M21_gp*migallow, #human population dynamics 
    D(P₂_gp) ~ r₂*P₂_gp*(1 - P₂_gp/Kp₂) - M21_gp, # standard logistic.
    D(K₁_gp) ~ s₁*Yd₁_gp - δ₁*K₁_gp - σ₁*exp(G_gp-G₁)*K₁_gp, #change in capital stock = 
    D(K₂_gp) ~ s₂*Yd₂_gp - δ₂*K₂_gp - σ₂*exp(G_gp-G₁)*K₂_gp, # savings-depreciation-damages
    D(G_gp)  ~ er₁_gp*Y₁_gp + er₂*Y₂ - u*(G_gp-2.8) + α*θ(G_gp-G₀), # global externality.
    D(z_gp)  ~ η*θ(G_gp-Tg)*(1-z_gp), #initiate decarbonization
    D(e₁_gp) ~ -z_gp*re₁*e₁_gp, # decarbonization in region 1
    D(e₂_gp) ~ -z_gp*re₂*e₂_gp, # decarbonization in region 2
    ## P1_gpc gdp pop and carrying capactiy
    D(P_kp) ~ r₂*P_kp*(1 - P_kp/kp_2),
    ## reference for carrying capacity driver 
    ########
    D(P₁_gpc) ~ r₁*P₁_gpc*(1 - P₁_gpc/Kp₁) + M21_gpc*migallow, #human population dynamics 
    D(P₂_gpc) ~ r₂*P₂_gpc*(1 - P₂_gpc/Kp₂) - M21_gpc, # standard logistic.
    D(K₁_gpc) ~ s₁*Yd₁_gpc - δ₁*K₁_gpc - σ₁*exp(G_gpc-G₁)*K₁_gpc, #change in capital stock = 
    D(K₂_gpc) ~ s₂*Yd₂_gpc - δ₂*K₂_gpc - σ₂*exp(G_gpc-G₁)*K₂_gpc, # savings-depreciation-damages
    D(G_gpc)  ~ er₁_gpc*Y₁_gpc + er₂_gpc*Y₂_gpc - u*(G_gpc-2.8) + α*θ(G_gpc-G₀), # global externality.
    D(z_gpc)  ~ η*θ(G_gpc-Tg)*(1-z_gpc), #initiate decarbonization
    D(e₁_gpc) ~ -z_gpc*re₁*e₁_gpc, # decarbonization in region 1
    D(e₂_gpc) ~ -z_gpc*re₂*e₂_gpc, # decarbonization in region 2
    ## all with remittances gpcr 
    D(P₁_gpcr) ~ r₁*P₁_gpcr*(1 - P₁_gpcr/Kp₁) + M21_gpcr*migallow, #human population dynamics 
    D(P₂_gpcr) ~ r₂*P₂_gpcr*(1 - P₂_gpcr/Kp₂) - M21_gpcr, # standard logistic.
    D(K₁_gpcr) ~ s₁*Yd₁_gpcr - δ₁*K₁_gpcr - σ₁*exp(G_gpcr-G₁)*K₁_gpcr - r*(D_21_gpcr*(Y₁_gpcr/P₁_gpcr)), #change in capital stock = 
    D(K₂_gpcr) ~ s₂*Yd₂_gpcr - δ₂*K₂_gpcr - σ₂*exp(G_gpcr-G₁)*K₂_gpcr + r*(D_21_gpcr*(Y₁_gpcr/P₁_gpcr)), # savings-depreciation-damages
    D(G_gpcr)  ~ er₁_gpcr*Y₁_gpcr + er₂_gpcr*Y₂_gpcr - u*(G_gpcr-2.8) + α*θ(G_gpcr-G₀), # global externality.
    D(z_gpcr)  ~ η*θ(G_gpcr-Tg)*(1-z_gpcr), #initiate decarbonization
    D(e₁_gpcr) ~ -z_gpcr*re₁*e₁_gpcr, # decarbonization in region 1
    D(e₂_gpcr) ~ -z_gpcr*re₂*e₂_gpcr, # decarbonization in region 2
    ### extending by variables without migration 
    D(P1_nomig) ~ r₁*P1_nomig*(1 - P1_nomig/Kp₁),
    D(P2_nomig) ~ r₂*P2_nomig*(1 - P2_nomig/Kp₂),
    D(K1_nomig) ~ s₁*Yd₁_nomig - δ₁*K1_nomig - σ₁*exp(G_nomig-G₁)*K1_nomig,
    D(K2_nomig) ~ s₂*Yd₂_nomig - δ₂*K2_nomig - σ₂*exp(G_nomig-G₁)*K2_nomig,
    D(G_nomig) ~ er₁_nomig*Y₁_nomig + er₂_nomig*Y₂_nomig - u*(G_nomig-2.8) + α*θ(G_nomig-G₀),
    D(z_nomig) ~ η*θ(G_nomig-Tg)*(1-z_nomig),
    D(e1_nomig) ~ -z_nomig*re₁*e1_nomig,
    D(e2_nomig)~ -z_nomig*re₂*e2_nomig,
    D(D12)~0,
    D(D21)~0
    ]

@named sys = ODESystem(eqs,t,vars,pars);



prob = ODEProblem(sys, [], (0.0,1000), [])
sol = solve(prob,saveat=1.0)

# x = y for equal distribution of income 
x = LinRange(1,10,10)
y = LinRange(1,10,10)

#per capita gni
plot(sol[Y₁_nomig/P1_nomig],color = "dodgerblue2", label = "LICs",xticks = 0:200:1000);
plot!(sol[Y₂_nomig/P2_nomig],color = "orangered", label = "HICs", ylims= (0,40));

plot!(sol[Y₁_gp/P₁_gp],color = "dodgerblue2", label = false,linestyle=:dashdotdot);
plot!(sol[Y₂_gp/P₂_gp],color = "orangered", label = false, ylims= (0,40),linestyle=:dashdotdot);

# plot!(sol[Y₁_gpc/P₁_gpc],color = "dodgerblue2", label = false,linestyle=:dashdot);
# plot!(sol[Y₂_gpc/P₂_gpc],color = "orangered", label = false, ylims= (0,40),linestyle=:dashdot);

# plot!(sol[Y₁_gpcr/P₁_gpcr],color = "dodgerblue2", label = false,linestyle=:dot);
# plot!(sol[Y₂_gpcr/P₂_gpcr],color = "orangered", label = false, ylims= (0,40),linestyle=:dot);

plot!(sol[Y₁/P₁],color = "dodgerblue2", label = false,linestyle=:dash);
plot!(sol[Y₂/P₂],xlabel="Time", ylabel="\$/10³",color = "orangered", label = false, linestyle=:dash,
title = "Per Capita GNI and Externality", legend_column = -1,legend_position = :topleft,  ylims= (0,40));

plot!(twinx(), sol,idxs=[100*G_nomig],color="palegreen4", label = "CO2",ylims = (0,850));
plot!(twinx(), sol,idxs=[100*G_gp],color="palegreen4", label = false,ylims = (0,850),linestyle=:dashdotdot);
# plot!(twinx(), sol,idxs=[100*G_gpc],color="palegreen4", label = false,ylims = (0,850),linestyle=:dashdot);
# plot!(twinx(), sol,idxs=[100*G_gpcr], color="palegreen4", label = false,ylims = (0,850),linestyle=:dot);

p1=plot!(twinx(), sol,idxs=[100*G],xlabel="", ylabel="Atm CO2(ppm)", color="palegreen4", label = false,
legend_position = :topright,  ylims = (0,850),linestyle=:dash);

#population

plot(sol[P1_nomig],color = "dodgerblue2", label = "HICs");
plot!(sol[P2_nomig],color = "orangered", label = "LICs");
plot!(sol[P2_nomig+P1_nomig],xlabel=false, color="grey", label = "World");

plot!(sol[P₁_gp],color = "dodgerblue2", label = false,linestyle=:dashdotdot);
plot!(sol[P₂_gp],color = "orangered", label = false,linestyle=:dashdotdot);
plot!(sol[P₁_gp+P₂_gp],color="grey", label = false,linestyle=:dashdotdot);

# plot!(sol[P₁_gpc],color = "dodgerblue2", label = false,linestyle=:dashdot);
# plot!(sol[P₂_gpc],color = "orangered", label = false, ylims= (0,40),linestyle=:dashdot);
# plot!(sol[P₂_gpc+P₁_gpc],color="grey", label = false,linestyle=:dashdot);

# plot!(sol[P₁_gpcr],color = "dodgerblue2", label = false,linestyle=:dot);
# plot!(sol[P₂_gpcr],color = "orangered", label = false, ylims= (0,40),linestyle=:dot);
# plot!(sol[P₂_gpcr+P₁_gpcr], color="grey", label = false,linestyle=:dot);

plot!(sol,idxs=[P₁],color = "dodgerblue2", label = false,linestyle=:dash);
plot!(sol,idxs=[P₂],color = "orangered", label = false,linestyle=:dash);
p2=plot!(sol[P₁+P₂],xlabel="Time", ylabel="Billion People", label=false, 
title = "Population",legend_column = -1,legend_position = :topleft, yticks = 0:5:20,
color = "gray", ylims= (0,15),linestyle=:dash); 

#Distribution 
plot(x,y,color = "orange", label = false);
plot!(sol[Y₁_nomig/P1_nomig],sol[Y₂_nomig/P2_nomig], color = "black", label = false);
plot!(sol[Y₁_gp/P₁_gp],sol[Y₂_gp/P₂_gp], color = "black", label = false,linestyle=:dashdotdot);
# plot!(sol[Y₁_gpc/P₁_gpc],sol[Y₂_gpc/P₂_gpc], color = "black", label = false,linestyle=:dashdot);
# plot!(sol[Y₁_gpcr/P₁_gpcr],sol[Y₂_gpcr/P₂_gpcr], color = "black", label = false,linestyle=:dot);
p3=plot!(sol[Y₁/P₁],sol[Y₂/P₂], xlabel="Per Capita GNI, HICs (\$/10³)",linestyle=:dash,
ylabel="Per Capita GNI, LICs (\$/10³)", color = "black", label = false, title = "JOS");


#operating space
plot(sol[G_nomig*100],sol[Y₁_nomig+Y₂_nomig], color= "black", label = false);
plot!(sol[G_gp*100],sol[Y₁_gp+Y₂_gp], color= "black", label = false,linestyle=:dashdotdot);
# plot!(sol[G_gpc*100],sol[Y₁_gpc+Y₂_gpc], color= "black", label = false,linestyle=:dashdotdot);
# plot!(sol[G_gpcr*100],sol[Y₁_gpcr+Y₂_gpcr], color= "black", label = false,linestyle=:dashdotdot);
p4=plot!(sol[G*100],sol[Y₁+Y₂], xlabel="Atm CO2(ppm)",
ylabel="Total World GNI(\$/10¹²)", color= "black", label = false, title="SOS",linestyle=:dash);

plt = plot(p2, p1, p3, p4, layout=(2,2), size=(750,500))
display(plt)
#safe

plot((sol[(Y₁*P₂)/(P₁*Y₂)]),label="a = 0.04, b = c = 0, rem = 0")
plot!((sol[(Y₁_gp*P₂_gp)/(P₁_gp*Y₂_gp)]),label="a = b 0.04, c = 0, rem = 0")
plot!((sol[(Y₁_gpc*P₂_gpc)/(P₁_gpc*Y₂_gpc)]),label="a = b = 0.04, c = 0.1, rem = 0")
#p5 = plot!((sol[Y₁_nomig/P1_nomig]/sol[Y₂_nomig/P2_nomig]), title = "PC1/PC2")
p5 = plot!((sol[(Y₁_nomig*P2_nomig)/ (P1_nomig*Y₂_nomig)]), title = "PC1/PC2", label = "nomig")

plot(sol[G*100],sol[P₁+P₂],label="a = 0.04, b=c = 0, rem = 0.1")
plot!(sol[G_gp*100],sol[P₁_gp+P₂_gp],label="a = b 0.04, c = 0, rem = 0.1")
plot!(sol[G_gpc*100],sol[P₁_gpc+P₂_gpc],label="a = b 0.04, c = 0.1, rem = 0.1")
p6 = plot!(sol[G_nomig]*100,sol[P1_nomig+P2_nomig], label = "nomig", title = "Externality vs Popuation", xlabel = "atm CO2 (ppm)", ylabel = "total world pop (billion ppl)")

plot(sol[(G*100)/(P₁+P₂)],label="a = 0.04, b=c = 0, rem = 0.1")
plot!(sol[(G_gp*100)/(P₁_gp+P₂_gp)],label="a = b 0.04, c = 0, rem = 0.1")
plot!(sol[(G_gpc*100)/(P₁_gpc+P₂_gpc)],label="a = b 0.04, c = 0.1, rem = 0.1")
p7 = plot!(sol[(G_nomig*100)/(P1_nomig+P2_nomig)], label = "nomig", title = "Externality/Popuation", xlabel = "time", ylabel = "G/worldpop")




plt2 = plot(p5,p6,p7 ,layout = (2,2))
display(plt2)



#plot(sol[D12])
#p5=plot!(sol[D12], xlabel=")",
#ylabel="", color= "black", title="Diaspora");

#savefig(plt,"test.png")

#colors: 1(HIC) color = "dodgerblue2", 2 (LIC)color = "orangered", CO2  color="palegreen4"