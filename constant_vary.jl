using ModelingToolkit, OrdinaryDiffEq, Plots, DifferentialEquations

#for i in 1:10

#m12 = LinRange(0)

@parameters t 
D = Differential(t)

parspop = @parameters r₁=0.038 Kp₁=1.5 r₂=0.042 Kp₂=9.7 mₒ=5 m₁=5 #population/migration dynamics
parsmig = @parameters r = 0.01 b1= 0.01 b2= 0.01 b3= 0.01 m12= 0.00 m21=0.004 #migration parameter
parsecn = @parameters δ₁=0.05 δ₂=0.05 α₁=0.5 α₂=0.5 s₁=0.25 s₂=0.21 Cₘ₁=0.7 Cₘ₂=0.7 a₁=2.7 a₂=1.7 #economic dynamics
parsear = @parameters r = 0.5 re₁=0.1 re₂=0.1 eb=0.00004 u=0.0025 α=0.1 σ₁=0.03 σ₂=0.03 G₀=20 G₁=5 Tg=50 η=1 #earth system and climate 
pars = [parspop;parsmig;parsecn;parsear] #combine into a single vector.

vars = @variables P₁(t)=0.24 P₂(t)=0.24 K₁(t)=0 K₂(t)=0 G(t)=2.8 e₁(t)=0.0004 e₂(t)=0.0004 z(t)=0 P1_nomig(t)=0.24 P2_nomig(t)=0.24 K1_nomig(t)=0 K2_nomig(t)=0 G_nomig(t)=2.8 e1_nomig(t)=0.0004 e2_nomig(t)=0.0004 z_nomig(t)=0


θ(x) = (max(0,x))^5/(0.0001^5 + (max(0,x))^5) #jump function 
@register θ(x)

ϕ(x) = x*θ(x) # linear for positive, zero for negative values 
@register ϕ(x)

Yᵢ(a,α,K,P) = a*((ϕ(K))^α)*((ϕ(P))^(1-α)) # Cobb-Douglas 
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

Mij(Pi,mij) = Pi.*mij
@register Mij(Pi,mij)

M12= Mij(P₁,m12)
M21= Mij(P₂,m21)

migallow= 0

D = Differential(t)
# r₁*P₁*p₁₂*(1 - P₁*p₁₂/Kp₁) + r₂*P₂*p₂₁*(1 - P₂*p₂₁/Kp₂)
#+ r₁*P₁*p₁₂*(1 - P₁*p₁₂/Kp₁) - r₂*P₂*p₂₁*(1 - P₂*p₂₁/Kp₂)
eqs = [
    D(P₁) ~ r₁*(P₁-M12)*(1 - (P₁-M12)/Kp₁) + M21,#+ r₂*M21*(1 - M21/Kp₂) *migallow, #human population dynamics 
    D(P₂) ~ r₂*(P₂-M21)*(1 - (P₂-M21)/Kp₂), #+ r₁*(M12)*(1-(M12)/Kp₁) , # standard logistic.
    D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - σ₁*exp(G-G₁)*K₁, #- r*(D21*(Y₁/P₁)), #change in capital stock = 
    D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - σ₂*exp(G-G₁)*K₂,# + r*(D21*(Y₁/P₁)), # savings-depreciation-damages
    D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8), #+ α*θ(G-G₀), # global externality.
    D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
    D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
    D(e₂) ~ -z*re₂*e₂, # decarbonization in region 2
    D(P1_nomig) ~ r₁*P1_nomig*(1 - P1_nomig/Kp₁),
    D(P2_nomig) ~ r₂*P2_nomig*(1 - P2_nomig/Kp₂),
    D(K1_nomig) ~ s₁*Yd₁_nomig - δ₁*K1_nomig - σ₁*exp(G_nomig-G₁)*K1_nomig,
    D(K2_nomig) ~ s₂*Yd₂_nomig - δ₂*K2_nomig - σ₂*exp(G_nomig-G₁)*K2_nomig,
    D(G_nomig) ~ er₁_nomig*Y₁_nomig + er₂_nomig*Y₂_nomig - u*(G_nomig-2.8), #+ α*θ(G_nomig-G₀),
    D(z_nomig) ~ η*θ(G_nomig-Tg)*(1-z_nomig),
    D(e1_nomig) ~ -z_nomig*re₁*e1_nomig,
    D(e2_nomig)~ -z_nomig*re₂*e2_nomig
    ]

@named sys = ODESystem(eqs,t,vars,pars);


prob = ODEProblem(sys, [], (0.0,1000), [])
sol = solve(prob,saveat=1.0)
    
# x = y for equal distribution of income 
x = LinRange(1,10,10)
y = LinRange(1,10,10)

P1 = sol[1,:]
P2 = sol[2,:]
k1 = sol[3,:]
k2 = sol[4,:]
YPC1 = k1./P1
YPC2 = k2./P2
GDP1 = 1 ./(1 .+YPC1./(YPC1+YPC2))
GDP2 = 1 ./(1 .+YPC2./(YPC2+YPC1))


    P1= sol[1,:];
    P2= sol[2,:];
    K1= sol[3,:];
    K2= sol[4,:];
    GG= sol[5,:];
    #em₁ = sol[7,:];
    #em₂ = sol[8,:];
    t = sol.t;

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

plt = plot(p2, p1, p3, p4, layout=(2,2), size=(750,500),legend =false)#legend_position=:outerbottom)
display(plt)



p5 = plot(sol[P₁] + sol[P₂],sol[G], xlabel="Total World Population",
ylabel="Atm CO2(ppm)", color = "black", title = "Carbon Population Ratio")
display(p5)


p6 = plot((sol[Y₁/P₁]./sol[Y₂/P₂]), xlabel="Time",
ylabel="Per Capita Ration LICs/HICs", color = "black", title = "Per Capita Income Ratio")
display(p6)

