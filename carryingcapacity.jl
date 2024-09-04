
using ModelingToolkit, DifferentialEquations, Plots


@parameters t
D = Differential(t)


parspop = @parameters r₁=0.038 Kp₁=1.5 r₂=0.042 Kp₂=9.7 mₒ=5 m₁=5 a = 0.1 b = 0.0005 τCO2 = 0.001 τPOP = 0.01#population/migration dynamics
parsecn = @parameters Ypc0 = 10 δ₁=0.05 δ₂=0.05 α₁=0.5 α₂=0.5 s₁=0.25 s₂=0.21 Cₘ₁=0.7 Cₘ₂=0.7 a₁=2.7 a₂=1.7 #economic dynamics
parsear = @parameters r = 0.5 re₁=0.1 re₂=0.1 eb=0.00004 u=0.0025 α=0.1 σ₁=0.03 σ₂=0.03 G₀=20 G₁=5 Tg=50 η=1 #earth system and climate 
pars = [parspop;parsecn;parsear] #combine into a single vector.

vars = @variables P₁(t)=0.24 P₂(t)=0.24 K₁(t)=0 K₂(t)=0 G(t)=2.8 e₁(t)=0.0004 e₂(t)=0.0004 z(t)=0 P1_nomig(t)=0.24 P2_nomig(t)=0.24 K1_nomig(t)=0 K2_nomig(t)=0 G_nomig(t)=2.8 e1_nomig(t)=0.0004 e2_nomig(t)=0.0004 z_nomig(t)=0 D12(t)=0 D21(t)=0  

θ(x) = (max(0,x))^5/(0.0001^5 + (max(0,x))^5)
@register θ(x)

ϕ(x) = x*θ(x)
@register ϕ(x)


Yᵢ(a,α,K,P) = a*((ϕ(K))^α)*((ϕ(P))^(1-α))
@register Yᵢ(a,α,K,P)

Yᵢ₁ = Yᵢ(a₁,α₁,K₁,P₁)
Yᵢ₂ = Yᵢ(a₂,α₂,K₂,P₂)
Yᵢ₁_nomig = Yᵢ(a₁,α₁,K1_nomig,P1_nomig)
Yᵢ₂_nomig = Yᵢ(a₂,α₂,K2_nomig,P2_nomig)

Y(P,Y) = P + ϕ(Y-P)
@register Y(P,K)

Y₁ = Y(P₁,Yᵢ₁)
Y₂ = Y(P₂,Yᵢ₂)
Y₁_nomig = Y(P1_nomig,Yᵢ₁_nomig)
Y₂_nomig = Y(P2_nomig,Yᵢ₂_nomig)

er₁ = θ(P₁-Yᵢ₁)*eb + θ(Yᵢ₁-P₁)*e₁
er₂ = θ(P₂-Yᵢ₂)*eb + θ(Yᵢ₂-P₂)*e₂
er₁_nomig = θ(P1_nomig-Yᵢ₁_nomig)*eb + θ(Yᵢ₁_nomig-P1_nomig)*e1_nomig
er₂_nomig = θ(P2_nomig-Yᵢ₂_nomig)*eb + θ(Yᵢ₂_nomig-P2_nomig)*e2_nomig

Yd₁ = ϕ(Y₁ - P₁*Cₘ₁)
Yd₂ = ϕ(Y₂ - P₂*Cₘ₂)
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
pop_i(Pi,Pj) = Pi./(Pi+Pj)
@register pop_i(Pi,Pj)

pop1 = pop_i(P₁,P₂) 
pop2 = pop_i(P₂,P₁)

#GDP difference function (share of per capita at home vs per capita destination)
gdp_ij(Yj,Pj,Yi,Pi) = (Yj./(Pj.*((Yj./Pj)+(Yi./Pi))))


## from more people to poor to move
gdp2_ij(Yj,Pj,Yi,Pi) = 1/(1+Yi./(Pi.*((Yj./Pj)+(Yi./Pi))))#.*(1/(1+exp.(γ*((Yi./Pi)-Ypc0))))
@register gdp2_ij(Yj,Pj,Yi,Pi)

#pij(Yj,Pj,Yi,Pi) = 1 .-(Yi./(Pi.*((Yi./Pi)+(Yj./Pj)))) 
#pij(Yj,Pj,Yi,Pi) = 1 -(Yj/(Pj*((Yj/Pj)+(Yi/Pi)))) ## falschrum
#@register pij(Yj,Pj,Yi,Pi)

kp(Kp,G) = Kp - G*θ(G-4)#(ϕ(G-4))
@register kp(Kp,G)
kp_2 = kp(Kp₂,G)

m12 = 0.0#gdp2_ij(Y₂,P₂,Y₁,P₁).*a - pop1.*b + 0.01*θ(P₁-KP1)
m21 = 0.01*θ(P₂-kp_2)#gdp2_ij(Y₁,P₁,Y₂,P₂).*a - pop2.*b + 0.01*θ(P₂-KP2)

Mij(P,m) = P.*m
@register Mij(P,m)

M12= Mij(P₁,m12)
M21= Mij(P₂,m21)

#define return migration different as push migration? 


#migallow, allows migrants to enetr or not 
#--> people disappear or die on the way or get stranded in between 
# --> system looses population
migallow = 1

#inconstant carrying capacity 

gamma(τC, τP, G,P,Kp)= τC*(G/280) + τP*(P/Kp)#τC*(G/280) + τP*(P/Kp)
@register gamma(τC, τP, G,P,Kp)
## difference in carrying capacity increases pressure to migrate -- when carrying capacity lowers --> habitability decreases

γ1 = gamma(0.001, 0.0,G,P₁,Kp₁)##add here the climate damage factor for both regions differently --> injust and different vulnerabiliy
γ2 = gamma(0.01, 0.01,G,P₂,Kp₂)


#differential equations

eqs = [
    D(P₁) ~ r₁*P₁*(1 - P₁/Kp₁) - M12 + M21*migallow, #human population dynamics 
    D(P₂) ~ r₂*P₂*(1 - P₂/Kp₂) + M12 - M21, # standard logistic.
    D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - σ₁*exp(G-G₁)*K₁ - r*(M12*(Y₁/P₁)), #change in capital stock = 
    D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - σ₂*exp(G-G₁)*K₂ + r*(M12*(Y₁/P₁)), # savings-depreciation-damages
    D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8) + α*θ(G-G₀), # global externality.
    D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
    D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
    D(e₂) ~ -z*re₂*e₂, # decarbonization in region 2
    ### extending by variables without migration 
    D(P1_nomig) ~ r₁*P1_nomig*(1 - P1_nomig/Kp₁),
    D(P2_nomig) ~ r₂*P2_nomig*(1 - P2_nomig/kp_2),
    D(K1_nomig) ~ s₁*Yd₁_nomig - δ₁*K1_nomig - σ₁*exp(G_nomig-G₁)*K1_nomig,
    D(K2_nomig) ~ s₂*Yd₂_nomig - δ₂*K2_nomig - σ₂*exp(G_nomig-G₁)*K2_nomig,
    D(G_nomig) ~ er₁_nomig*Y₁_nomig + er₂_nomig*Y₂_nomig - u*(G_nomig-2.8) + α*θ(G_nomig-G₀),
    D(z_nomig) ~ η*θ(G_nomig-Tg)*(1-z_nomig),
    D(e1_nomig) ~ -z_nomig*re₁*e1_nomig,
    D(e2_nomig)~ -z_nomig*re₂*e2_nomig,
    D(D12)~M12,
    D(D21)~M21,
    ]

@named sys = ODESystem(eqs,t,vars,pars);



prob = ODEProblem(sys, [], (0.0,1000), [])
sol = solve(prob,saveat=1.0)


# x = y for equal distribution of income 
x = LinRange(1,10,10)
y = LinRange(1,10,10)

#plot no migration

#population
plot(sol,idxs=[P1_nomig],color = "dodgerblue2", label = "HICs");
plot!(sol,idxs=[P2_nomig],color = "orangered", label = "LICs");
p2=plot!(sol[P1_nomig+P2_nomig],xlabel="Time", ylabel="billion people)", label="World", 
title = "Population",legend_column = -1,legend_position = :topleft, yticks = 0:5:20,
ylims= (0,15)); 

display(p2)

# #per capita gni
# plot(sol[Y₁/P₁],color = "dodgerblue2", label = "HICs");
# plot!(sol[Y₂/P₂],xlabel="Time", ylabel="\$/10³", 
# color = "orangered", label = "LICs", title = "Per capita GNI and Externality",
# legend_column = -1,legend_position = :topleft, yticks = 0:5:40, ylims= (0,40));
# p1=plot!(twinx(), sol,idxs=[10*G],xlabel="", 
# ylabel="Atmoshperic CO2(ppm)", color="palegreen4", label = "CO2",
# legend_position = :topright, yticks = 0:10:75, ylims = (0,85));

# #population
# plot(sol,idxs=[P₁],color = "dodgerblue2", label = "HICs");
# plot!(sol,idxs=[P₂],color = "orangered", label = "LICs");
# p2=plot!(sol[P₁+P₂],xlabel="Time", ylabel="billion people)", label="World", 
# title = "Population",legend_column = -1,legend_position = :topleft, yticks = 0:5:20,
# ylims= (0,15)); 
# plot(x,y,color = "orange", label = false);
# #Distribution 
# p3=plot!(sol[Y₁/P₁],sol[Y₂/P₂], xlabel="per capita GNI, HICs (\$/10³)",
# ylabel="per capita GNI, LICs (\$/10³)", color = "black", label = false, title = "GNI distribution");
# #operating space
# p4=plot(sol[G*100],sol[Y₁+Y₂], xlabel="Atmospheric CO2(ppm)",
# ylabel="Total World GNI(\$/10¹²)", color= "black", label = false, title="Operating Space");
# plt = plot(p2, p1, p3, p4, layout=(2,2), size=(750,500))

#display(plt)
#safe

# plot(sol[D12])
# p5=plot!(sol[D12], xlabel=")",
# ylabel="", color= "black", title="Diaspora");

#savefig(plt,"test.png")

#colors: 1(HIC) color = "dodgerblue2", 2 (LIC)color = "orangered", CO2  color="palegreen4"



###GNI per capita and Atmospheric CO2
plot(sol[Y₁/P₁],color = "dodgerblue2", label = "HICs");
plot!(sol[Y₁_nomig/P1_nomig],color = "dodgerblue2",linestyle=:dash,label=false);
plot!(sol[Y₂_nomig/P2_nomig],color = "orangered",linestyle=:dash,label=false);
plot!(sol[Y₂/P₂],xlabel="Time", ylabel="\$/10³",color = "orangered", label = "LICs", 
title = "Per capita GNI and Externality",yticks = 0:5:40, ylims= (0,40));

plot!(twinx(), sol,idxs=[10*G_nomig],xlabel="",color="palegreen4",linestyle=:dash,legend=false
,yticks = 0:10:75, ylims = (0,85))

p1=plot!(twinx(), sol,idxs=[10*G],xlabel="", ylabel="Atmospheric CO2(ppm)", color="palegreen4",
label = "CO2",yticks = 0:10:75, ylims = (0,85));

#population
plot(sol,idxs=[P1_nomig],color = "dodgerblue2",linestyle=:dash);
plot!(sol,idxs=[P2_nomig],color = "orangered",linestyle=:dash);
plot!(sol[P1_nomig+P2_nomig],xlabel="Time", ylabel="billion people)",
linestyle=:dash, color = "gray");
plot!(sol,idxs=[P₁],color = "dodgerblue2", label = "HICs");
plot!(sol,idxs=[P₂],color = "orangered", label = "LICs");
p2=plot!(sol[P₁+P₂],xlabel="Time", ylabel="billion people)", label="World", 
title = "Population",color="gray",yticks = 0:5:20,ylims= (0,15)); 

#Distribution
plot(x,y,color = "orangered", label = "equal distribution"); 
plot!(sol[Y₁_nomig/P1_nomig],sol[Y₂_nomig/P2_nomig], xlabel="per capita GNI, HICs (\$/10³)",
linestyle=:dash,color = "black");
p3=plot!(sol[Y₁/P₁],sol[Y₂/P₂], xlabel="per capita GNI, HICs (\$/10³)",
ylabel="per capita GNI, LICs (\$/10³)", color = "black", title = "GNI distribution")


#operating space
plot(sol[G_nomig*100],sol[Y₁_nomig+Y₂_nomig],linestyle=:dash,color="green");
p4=plot!(sol[G*100],sol[Y₁+Y₂], xlabel="Atmospheric CO2(ppm)",
ylabel="Total World GNI(\$/10¹²)", color= "green", title="Operating Space");

## carrying capacity 
#plot(sol,idxs=[KP1],color = "dodgerblue2",linestyle=:dash);
#p5 = plot!(sol,idxs=[KP2],color = "orangered",linestyle=:dash);

plt = plot(p2, p1, p3, p4, layout=(2,2), size=(750,500),legend =false)#legend_position=:outerbottom)
display(plt)


#savefig(plt,"carryingcapacity_theta.svg")