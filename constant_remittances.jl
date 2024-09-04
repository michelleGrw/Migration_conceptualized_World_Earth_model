using ModelingToolkit, OrdinaryDiffEq, Plots, DifferentialEquations


@parameters t
D = Differential(t)

parspop = @parameters r₁=0.038 Kp₁=1.5 r₂=0.042 Kp₂=9.7 mₒ=5 m₁=5 #population/migration dynamics
parsmig = @parameters r = 0.4 b1= 0.01 b2= 0.01 b3= 0.01 m12= 0.00 m21=0.01 #migration parameter
parsecn = @parameters δ₁=0.05 δ₂=0.05 α₁=0.5 α₂=0.5 s₁=0.25 s₂=0.21 Cₘ₁=0.7 Cₘ₂=0.7 a₁=2.7 a₂=1.7 #economic dynamics
parsear = @parameters r = 0.5 re₁=0.1 re₂=0.1 eb=0.00004 u=0.0025 α=0.1 σ₁=0.03 σ₂=0.03 G₀=20 G₁=5 Tg=50 η=1 #earth system and climate 
pars = [parspop;parsmig;parsecn;parsear] #combine into a single vector.

vars = @variables P₁(t)=0.24 P₂(t)=0.24 K₁(t)=0 K₂(t)=0 G(t)=2.8 e₁(t)=0.0004 e₂(t)=0.0004 z(t)=0 P₁_rem(t)=0.24 P₂_rem(t)=0.24 K₁_rem(t)=0 K₂_rem(t)=0 G_rem(t)=2.8 e₁_rem(t)=0.0004 e₂_rem(t)=0.0004 z_rem(t)=0 P₁_rem3(t)=0.24 P₂_rem3(t)=0.24 K₁_rem3(t)=0 K₂_rem3(t)=0 G_rem3(t)=2.8 e₁_rem3(t)=0.0004 e₂_rem3(t)=0.0004 z_rem3(t)=0 P₁_rem2(t)=0.24 P₂_rem2(t)=0.24 K₁_rem2(t)=0 K₂_rem2(t)=0 G_rem2(t)=2.8 e₁_rem2(t)=0.0004 e₂_rem2(t)=0.0004 z_rem2(t)=0 D12(t) = 0 D21(t) = 0 P1_nomig(t)=0.24 P2_nomig(t)=0.24 K1_nomig(t)=0 K2_nomig(t)=0 G_nomig(t)=2.8 e1_nomig(t)=0.0004 e2_nomig(t)=0.0004 z_nomig(t)=0 


θ(x) = (max(0,x))^5/(0.0001^5 + (max(0,x))^5) #jump function 
@register θ(x)

ϕ(x) = x*θ(x) # linear for positive, zero for negative values 
@register ϕ(x)

Yᵢ(a,α,K,P) = a*((ϕ(K))^α)*((ϕ(P))^(1-α)) # Cobb-Douglas 
@register Yᵢ(a,α,K,P)

Yᵢ₁ = Yᵢ(a₁,α₁,K₁,P₁)
Yᵢ₂ = Yᵢ(a₂,α₂,K₂,P₂)
Yᵢ₁_rem = Yᵢ(a₁,α₁,K₁_rem,P₁_rem)
Yᵢ₂_rem = Yᵢ(a₂,α₂,K₂_rem,P₂_rem)
Yᵢ₁_rem2 = Yᵢ(a₁,α₁,K₁_rem2,P₁_rem2)
Yᵢ₂_rem2 = Yᵢ(a₂,α₂,K₂_rem2,P₂_rem2)
Yᵢ₁_rem3 = Yᵢ(a₁,α₁,K₁_rem3,P₁_rem3)
Yᵢ₂_rem3 = Yᵢ(a₂,α₂,K₂_rem3,P₂_rem3)
Yᵢ₁_nomig = Yᵢ(a₁,α₁,K1_nomig,P1_nomig)
Yᵢ₂_nomig = Yᵢ(a₂,α₂,K2_nomig,P2_nomig)

Y(P,Y) = P + ϕ(Y-P)
@register Y(P,K)

Y₁ = Y(P₁,Yᵢ₁)
Y₂ = Y(P₂,Yᵢ₂)
Y₁_rem = Y(P₁_rem,Yᵢ₁_rem)
Y₂_rem = Y(P₂_rem,Yᵢ₂_rem)
Y₁_rem2 = Y(P₁_rem2,Yᵢ₁_rem2)
Y₂_rem2 = Y(P₂_rem2,Yᵢ₂_rem2)
Y₁_rem3 = Y(P₁_rem3,Yᵢ₁_rem3)
Y₂_rem3 = Y(P₂_rem3,Yᵢ₂_rem3)
Y₁_nomig = Y(P1_nomig,Yᵢ₁_nomig)
Y₂_nomig = Y(P2_nomig,Yᵢ₂_nomig)

er₁ = θ(P₁-Yᵢ₁)*eb + θ(Yᵢ₁-P₁)*e₁
er₂ = θ(P₂-Yᵢ₂)*eb + θ(Yᵢ₂-P₂)*e₂
er₁_nomig = θ(P1_nomig-Yᵢ₁_nomig)*eb + θ(Yᵢ₁_nomig-P1_nomig)*e1_nomig
er₂_nomig = θ(P2_nomig-Yᵢ₂_nomig)*eb + θ(Yᵢ₂_nomig-P2_nomig)*e2_nomig
er₁_rem = θ(P₁_rem-Yᵢ₁_rem)*eb + θ(Yᵢ₁_rem-P₁_rem)*e₁_rem
er₂_rem = θ(P₂_rem-Yᵢ₂_rem)*eb + θ(Yᵢ₂_rem-P₂_rem)*e₂_rem
er₁_rem2 = θ(P₁_rem2-Yᵢ₁_rem2)*eb + θ(Yᵢ₁_rem2-P₁_rem2)*e₁_rem2
er₂_rem2 = θ(P₂_rem2-Yᵢ₂_rem2)*eb + θ(Yᵢ₂_rem2-P₂_rem2)*e₂_rem2
er₁_rem3 = θ(P₁_rem3-Yᵢ₁_rem3)*eb + θ(Yᵢ₁_rem3-P₁_rem3)*e₁_rem3
er₂_rem3 = θ(P₂_rem3-Yᵢ₂_rem3)*eb + θ(Yᵢ₂_rem3-P₂_rem3)*e₂_rem3

Yd₁ = ϕ(Y₁ - P₁*Cₘ₁)
Yd₂ = ϕ(Y₂ - P₂*Cₘ₂)
Yd₁_rem = ϕ(Y₁_rem - P₁_rem*Cₘ₁)
Yd₂_rem = ϕ(Y₂_rem - P₂_rem*Cₘ₂)
Yd₁_rem2 = ϕ(Y₁_rem2 - P₁_rem2*Cₘ₁)
Yd₂_rem2 = ϕ(Y₂_rem2 - P₂_rem2*Cₘ₂)
Yd₁_rem3 = ϕ(Y₁_rem3 - P₁_rem3*Cₘ₁)
Yd₂_rem3 = ϕ(Y₂_rem3 - P₂_rem3*Cₘ₂)
Yd₁_nomig = ϕ(Y₁_nomig - P1_nomig*Cₘ₁)
Yd₂_nomig = ϕ(Y₂_nomig - P2_nomig*Cₘ₂)

Mij(Pi,mij) = Pi.*mij
@register Mij(Pi,mij)

M12= Mij(P₁,m12)
M21= Mij(P₂,m21)
M21_rem=Mij(P₂_rem,m21)
M21_rem2=Mij(P₂_rem2,m21)
M21_rem3=Mij(P₂_rem3,m21)

Dij(Mij,ri,Kpi) = Mij*ri*(1 - (Mij)/Kpi)
@register Dij(Mij,ri,Kpi)

D21_ind = P₁_rem-P1_nomig#Dij(M21_rem,r₂,Kp₂)
D21_ind2 = P₁_rem2-P1_nomig#Dij(M21_rem,r₂,Kp₂)
D21_ind3 = P₁_rem3-P1_nomig#Dij(M21_rem,r₂,Kp₂)

migallow= 0

rem = LinRange(0.05,0.2,4)

D = Differential(t)
# r₁*P₁*p₁₂*(1 - P₁*p₁₂/Kp₁) + r₂*P₂*p₂₁*(1 - P₂*p₂₁/Kp₂)
#+ r₁*P₁*p₁₂*(1 - P₁*p₁₂/Kp₁) - r₂*P₂*p₂₁*(1 - P₂*p₂₁/Kp₂)
eqs = [
    D(P₁) ~ r₁*(P₁)*(1 - (P₁/Kp₁)) + M21,#human population dynamics 
    D(P₂) ~ r₂*(P₂-M21)*(1 - (P₂-M21)/Kp₂), # standard logistic.
    D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - σ₁*exp(G-G₁)*K₁, #change in capital stock = 
    D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - σ₂*exp(G-G₁)*K₂,# savings-depreciation-damages
    D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8) + α*θ(G-G₀), # global externality.
    D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
    D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
    D(e₂) ~ -z*re₂*e₂,
    D(P₁_rem) ~ r₁*(P₁_rem)*(1 - (P₁_rem/Kp₁)) + M21_rem,#+ r₂*M21*(1 - M21/Kp₂) *migallow, #human population dynamics 
    D(P₂_rem) ~ r₂*(P₂_rem-M21_rem)*(1 - (P₂_rem-M21_rem)/Kp₂), #+ r₁*(M12)*(1-(M12)/Kp₁) , # standard logistic.
    D(K₁_rem) ~ s₁*Yd₁_rem - δ₁*K₁_rem - σ₁*exp(G_rem-G₁)*K₁_rem - rem[1]*(D21_ind*(Y₁_rem/P₁_rem)), #change in capital stock = 
    D(K₂_rem) ~ s₂*Yd₂_rem - δ₂*K₂_rem - σ₂*exp(G_rem-G₁)*K₂_rem + rem[1]*(D21_ind*(Y₁_rem/P₁_rem)), # savings-depreciation-damages
    D(G_rem)  ~ er₁*Y₁_rem + er₂*Y₂_rem - u*(G_rem-2.8) + α*θ(G_rem-G₀), # global externality.
    D(z_rem)  ~ η*θ(G_rem-Tg)*(1-z_rem), #initiate decarbonization
    D(e₁_rem) ~ -z_rem*re₁*e₁_rem, # decarbonization in region 1
    D(e₂_rem) ~ -z_rem*re₂*e₂_rem, # decarbonization in region 2
    D(P₁_rem2) ~ r₁*(P₁_rem)*(1 - (P₁_rem/Kp₁)) + M21_rem,#+ r₂*M21*(1 - M21/Kp₂) *migallow, #human population dynamics 
    D(P₂_rem2) ~ r₂*(P₂_rem-M21_rem)*(1 - (P₂_rem-M21_rem)/Kp₂), #+ r₁*(M12)*(1-(M12)/Kp₁) , # standard logistic.
    D(K₁_rem2) ~ s₁*Yd₁_rem2 - δ₁*K₁_rem2 - σ₁*exp(G_rem2-G₁)*K₁_rem2 - rem[2]*(D21_ind2*(Y₁_rem/P₁_rem)), #change in capital stock = 
    D(K₂_rem2) ~ s₂*Yd₂_rem2 - δ₂*K₂_rem2 - σ₂*exp(G_rem2-G₁)*K₂_rem2 + rem[2]*(D21_ind2*(Y₁_rem/P₁_rem)), # savings-depreciation-damages
    D(G_rem2)  ~ er₁*Y₁_rem2 + er₂*Y₂_rem2 - u*(G_rem2-2.8) + α*θ(G_rem-G₀), # global externality.
    D(z_rem2)  ~ η*θ(G_rem2-Tg)*(1-z_rem2), #initiate decarbonization
    D(e₁_rem2) ~ -z_rem2*re₁*e₁_rem2, # decarbonization in region 1
    D(e₂_rem2) ~ -z_rem2*re₂*e₂_rem2, # decarbonization in region 2
    D(P₁_rem3) ~ r₁*(P₁_rem3)*(1 - (P₁_rem3/Kp₁)) + M21_rem3,#+ r₂*M21*(1 - M21/Kp₂) *migallow, #human population dynamics 
    D(P₂_rem3) ~ r₂*(P₂_rem3-M21_rem3)*(1 - (P₂_rem3-M21_rem3)/Kp₂), #+ r₁*(M12)*(1-(M12)/Kp₁) , # standard logistic.
    D(K₁_rem3) ~ s₁*Yd₁_rem3 - δ₁*K₁_rem3 - σ₁*exp(G_rem3-G₁)*K₁_rem3 - rem[4]*(D21_ind3*(Y₁_rem3/P₁_rem3)), #change in capital stock = 
    D(K₂_rem3) ~ s₂*Yd₂_rem3 - δ₂*K₂_rem3 - σ₂*exp(G_rem3-G₁)*K₂_rem3 + rem[4]*(D21_ind3*(Y₁_rem3/P₁_rem3)), # savings-depreciation-damages
    D(G_rem3)  ~ er₁*Y₁_rem3 + er₂*Y₂_rem3 - u*(G_rem3-2.8) + α*θ(G_rem3-G₀), # global externality.
    D(z_rem3)  ~ η*θ(G_rem3-Tg)*(1-z_rem3), #initiate decarbonization
    D(e₁_rem3) ~ -z_rem3*re₁*e₁_rem3, # decarbonization in region 1
    D(e₂_rem3) ~ -z_rem3*re₂*e₂_rem3, # decarbonization in region 2
    D(D12) ~ 0, # diaspora only migration and return (migration from i only assumed from "natives")
    D(D21) ~ M21_rem*r₂*(1 - (M21_rem)/Kp₂), ### extending by variables without migration 
    D(P1_nomig) ~ r₁*P1_nomig*(1 - P1_nomig/Kp₁),
    D(P2_nomig) ~ r₂*P2_nomig*(1 - P2_nomig/Kp₂),
    D(K1_nomig) ~ s₁*Yd₁_nomig - δ₁*K1_nomig - σ₁*exp(G_nomig-G₁)*K1_nomig,
    D(K2_nomig) ~ s₂*Yd₂_nomig - δ₂*K2_nomig - σ₂*exp(G_nomig-G₁)*K2_nomig,
    D(G_nomig) ~ er₁_nomig*Y₁_nomig + er₂_nomig*Y₂_nomig - u*(G_nomig-2.8)+ α*θ(G_nomig-G₀),
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

### pc income
    plot(sol[Y₁/P₁],color = "dodgerblue2",linestyle=:dot, label = false);
    plot!(sol[Y₁_nomig/P1_nomig],color = "dodgerblue2",label="HICs");
    plot!(sol[Y₁_rem/P₁_rem],color = "dodgerblue2",linestyle=:dash,label=false);
    plot!(sol[Y₁_rem2/P₁_rem2],color = "dodgerblue2",linestyle=:dash,label=false);
    plot!(sol[Y₁_rem3/P₁_rem3],color = "dodgerblue2",linestyle=:dash,label=false);

    plot!(sol[Y₂_rem/P₂_rem],color = "orangered",label=false,linestyle=:dash);
    plot!(sol[Y₂_rem2/P₂_rem2],color = "orangered",label=false,linestyle=:dash);
    plot!(sol[Y₂_rem3/P₂_rem3],color = "orangered",label=false,linestyle=:dash);
    plot!(sol[Y₂_nomig/P2_nomig],color = "orangered",label="LICs");
    plot!(sol[Y₂/P₂],xlabel="Time", ylabel="\$/10³",color = "orangered", label = false, 
    title = "Per capita GNI and Externality",yticks = 0:10:30, ylims= (0,30),linestyle=:dot);
    
    plot!(twinx(), sol[100*G_nomig],xlabel="",color="green",legend=false,ylims = (0,700))
    plot!(twinx(), sol[100*G_rem],xlabel="",color="green",linestyle=:dash,legend=false,ylims = (0,700))
    plot!(twinx(), sol[100*G_rem2],xlabel="",color="green",linestyle=:dash,legend=false,ylims = (0,700))
    plot!(twinx(), sol[100*G_rem3],xlabel="",color="green",linestyle=:dash,legend=false,ylims = (0,700))
    p1=plot!(twinx(), sol,idxs=[100*G],xlabel="", ylabel="Atm CO2(ppm)", color="green",ylims = (0,700),
    label = "CO2", linestyle=:dot);
    
    # #population
    plot(sol,idxs=[P1_nomig],color = "dodgerblue2",label="HICs")#,leg=:outertop,legend_column = -1);
    plot!(sol,idxs=[P2_nomig],color = "orangered",label="LICs");
    plot!(sol[P1_nomig+P2_nomig],color = "grey",label="World");
    #plot!(legend=:outerbottom, legendcolumns=3)

    plot!(sol,idxs=[P₁],color = "dodgerblue2",linestyle=:dot, label = false);
    plot!(sol,idxs=[P₂],color = "orangered",linestyle=:dot, label = false);
    plot!(sol[P₁+P₂],color = "grey",linestyle=:dot, label = false);

    plot!(sol,idxs=[P₁_rem],color = "dodgerblue2", linestyle=:dash, label = false);
    plot!(sol,idxs=[P₂_rem],color = "orangered", linestyle=:dash, label = false);
    plot!(sol[P₁_rem+P₂_rem],xlabel="Time", ylabel="Billion People",
    title = "Population",yticks = 0:5:20,ylims= (0,15),color = "grey",linestyle=:dash, label = false); 

    plot!(sol,idxs=[P₁_rem2],color = "dodgerblue2",linestyle=:dash, label = false);
    plot!(sol,idxs=[P₂_rem2],color = "orangered", linestyle=:dash, label = false);
    plot!(sol[P₁_rem2+P₂_rem2],xlabel="Time", ylabel="Billione People",
    title = "Population",yticks = 0:5:20,ylims= (0,15),color = "grey",linestyle=:dash, label = false); 

    plot!(sol,idxs=[P₁_rem3],color = "dodgerblue2", linestyle=:dash, label = false);
    plot!(sol,idxs=[P₂_rem3],color = "orangered",linestyle=:dash, label = false);
    p2=plot!(sol[P₁_rem3+P₂_rem3],xlabel="Time", ylabel="Billione People",
    title = "Population",yticks = 0:5:20,ylims= (0,15),color = "grey",linestyle=:dash, label = false); 
    
    #JOS

    plot(x,y,color = "orangered", label = false); 
    plot!(sol[Y₁_nomig/P1_nomig],sol[Y₂_nomig/P2_nomig],color = "black", label = false);
    plot!(sol[Y₁_rem/P₁_rem],sol[Y₂_rem/P₂_rem],linestyle=:dash,color = "gray80", label = false);
    plot!(sol[Y₁_rem2/P₁_rem2],sol[Y₂_rem2/P₂_rem2],linestyle=:dash,color = "gray60", label = false);
    plot!(sol[Y₁_rem3/P₁_rem3],sol[Y₂_rem3/P₂_rem3],linestyle=:dash,color = "gray40", label = false);
    p3=plot!(sol[Y₁/P₁],sol[Y₂/P₂], xlabel="Per Capita GNI, HICs (\$/10³)",linestyle=:dot, label = false,
    ylabel="Per Capita GNI, LICs (\$/10³)", color = "black", title = "JOS")
    
    
    #SOS
    plot(sol[G_nomig*100],sol[Y₁_nomig+Y₂_nomig],color="black", label = false);
    plot!(sol[G_rem*100],sol[Y₁_rem+Y₂_rem],color="gray80",linestyle=:dash, label = false);
    plot!(sol[G_rem2*100],sol[Y₁_rem2+Y₂_rem2],linestyle=:dash,color="gray60", label = false);
    plot!(sol[G_rem3*100],sol[Y₁_rem3+Y₂_rem3],linestyle=:dash,color="gray40", label = false);
    p4=plot!(sol[G*100],sol[Y₁+Y₂], xlabel="Atm CO2(ppm)",linestyle=:dot,
    ylabel="Total World GNI(\$/10¹²)", color= "black", title="SOS", label = false);
    
    #p5 = plot(sol[D21], title = "Diaspora 2 in 1")
    plt = plot(p2, p1, p3, p4,layout=(2,2), size=(750,500),legend =false)

    display(plt)
    
    #savefig(plt,"remittances.svg")



# 