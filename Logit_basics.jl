
using ModelingToolkit, DifferentialEquations, Plots


@parameters t
D = Differential(t)


parspop = @parameters r₁=0.038 Kp₁=1.5 r₂=0.042 Kp₂=9.7 mₒ=5 m₁=5 #population/migration dynamics
parsmig = @parameters r = 0 b1= 0.02 b2= 0.02 b3= 0.0  a = 0#migration parameter
parsecn = @parameters δ₁=0.05 δ₂=0.05 α₁=0.5 α₂=0.5 s₁=0.25 s₂=0.21 Cₘ₁=0.7 Cₘ₂=0.7 a₁=2.7 a₂=1.7 #economic dynamics
parsear = @parameters re₁=0.1 re₂=0.1 eb=0.00004 u=0.0025 α=0.1 σ₁=0.03 σ₂=0.03 G₀=20 G₁=5 Tg=50 η=1 #earth system and climate 
pars = [parspop;parsmig;parsecn;parsear] #combine into a single vector.

vars = @variables P₁(t)=0.24 P₂(t)=0.24 K₁(t)=0 K₂(t)=0 G(t)=2.8 e₁(t)=0.0004 e₂(t)=0.0004 z(t)=0 P₁_gpc(t)=0.24 P₂_gpc(t)=0.24 K₁_gpc(t)=0 K₂_gpc(t)=0 G_gpc(t)=2.8 e₁_gpc(t)=0.0004 e₂_gpc(t)=0.0004 z_gpc(t)=0 P_kp(t)=0.24 P_kp1(t)=0.24 P1_nomig(t)=0.24 P2_nomig(t)=0.24 K1_nomig(t)=0 K2_nomig(t)=0 G_nomig(t)=2.8 e1_nomig(t)=0.0004 e2_nomig(t)=0.0004 z_nomig(t)=0

θ(x) = (max(0,x))^5/(0.0001^5 + (max(0,x))^5)
@register θ(x)

ϕ(x) = x*θ(x)
@register ϕ(x)

Yᵢ(a,α,K,P) = a*((ϕ(K))^α)*((ϕ(P))^(1-α))
@register Yᵢ(a,α,K,P)

Yᵢ₁ = Yᵢ(a₁,α₁,K₁,P₁)
Yᵢ₂ = Yᵢ(a₂,α₂,K₂,P₂)
Yᵢ₁_gpc = Yᵢ(a₁,α₁,K₁_gpc,P₁_gpc)
Yᵢ₂_gpc = Yᵢ(a₂,α₂,K₂_gpc,P₂_gpc)
Yᵢ₁_nomig = Yᵢ(a₁,α₁,K1_nomig,P1_nomig)
Yᵢ₂_nomig = Yᵢ(a₂,α₂,K2_nomig,P2_nomig)

Y(P,Y) = P + ϕ(Y-P)
@register Y(P,K)

Y₁ = Y(P₁,Yᵢ₁)
Y₂ = Y(P₂,Yᵢ₂)
Y₁_gpc = Y(P₁_gpc,Yᵢ₁_gpc)
Y₂_gpc = Y(P₂_gpc,Yᵢ₂_gpc)
Y₁_nomig = Y(P1_nomig,Yᵢ₁_nomig)
Y₂_nomig = Y(P2_nomig,Yᵢ₂_nomig)

er₁ = θ(P₁-Yᵢ₁)*eb + θ(Yᵢ₁-P₁)*e₁
er₂ = θ(P₂-Yᵢ₂)*eb + θ(Yᵢ₂-P₂)*e₂
er₁_gpc = θ(P₁_gpc-Yᵢ₁_gpc)*eb + θ(Yᵢ₁_gpc-P₁_gpc)*e₁_gpc
er₂_gpc = θ(P₂_gpc-Yᵢ₂_gpc)*eb + θ(Yᵢ₂_gpc-P₂_gpc)*e₂_gpc
er₁_nomig = θ(P1_nomig-Yᵢ₁_nomig)*eb + θ(Yᵢ₁_nomig-P1_nomig)*e1_nomig
er₂_nomig = θ(P2_nomig-Yᵢ₂_nomig)*eb + θ(Yᵢ₂_nomig-P2_nomig)*e2_nomig



#Diaspora(P,M) for remittances 
dia(P,P_nm)=P-P_nm
@register dia(P,P_nm)
D21 = dia(P₁,P1_nomig)
D21_gpc = dia(P₁_gpc,P1_nomig)

Yd₁ = ϕ(Y₁ - P₁*Cₘ₁ - r*(D21*(Y₁/P₁)))
Yd₂ = ϕ(Y₂ - P₂*Cₘ₂ + r*(D21*(Y₁/P₁)))
Yd₁_gpc = ϕ(Y₁_gpc - P₁_gpc*Cₘ₁) - r*(D21_gpc*(Y₁_gpc/P₁_gpc))
Yd₂_gpc = ϕ(Y₂_gpc - P₂_gpc*Cₘ₂) + r*(D21_gpc*(Y₁_gpc/P₁_gpc))
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

expo(b1,x1,b2,x2,b3,x3) = -5+b1*x1+b2*x2+b3*x3
@register expo(b1,x1,b2,x2,b3,x3)

logit(expo) = 1/(1+exp(-expo))
@register logit(expo)

pop_i(Pi,Pj) = (Pi./(Pj+Pi))
@register pop_i(Pi,Pj)

pop1 = pop_i(P₁,P₂) 
pop2 = pop_i(P₂,P₁)
pop2_gpc = pop_i(P₂_gpc,P₁_gpc)

#GDP difference function (share of per capita at home vs per capita destination)
ypc1 = Y₁./P₁
ypc2 = Y₂./P₂
ypc1_gpc = Y₁_gpc./P₁_gpc
ypc2_gpc = Y₂_gpc./P₂_gpc

## from more people to poor to move
gdp_ij(ypc1,ypc2) = (ypc1/(ypc2+ypc1))#.*(1/(1+exp.(γ*((Yi./Pi)-Ypc0)))) ##(ypc1/(ypc2+ypc1))
@register gdp_ij(ypc1,ypc2)

gdp1 = gdp_ij(ypc1,ypc2)
gdp2 = gdp_ij(ypc2,ypc1)
gdp2_gpc = gdp_ij(ypc2_gpc,ypc1_gpc)

kp(Kp,G) = Kp - G*θ(G-4)#(ϕ(G-4))
@register kp(Kp,G)
kp_1 = kp(Kp₂,G)
kp_2 = kp(Kp₂,G_gpc)


####this is carrying capacity factor that for expo function
kp_gp = θ(P₂-P_kp1)
kp_gpc = θ(P₂_gpc-P_kp)

expo21 = expo(5,gdp2,-3,pop2,0.5,kp_gp)
expo21_gpc = expo(3,gdp2_gpc,-2,pop2_gpc,0.5,kp_gpc)
expo12 = 0 #expo(1,0,-1,0,0.1,D21)

#pij(Yj,Pj,Yi,Pi) = 1 .-(Yi./(Pi.*((Yi./Pi)+(Yj./Pj)))) 
#pij(Yj,Pj,Yi,Pi) = 1 -(Yj/(Pj*((Yj/Pj)+(Yi/Pi)))) ## falschrum
#@register pij(Yj,Pj,Yi,Pi)

#logit(expo12)#logit(b1,gdp1,b2,pop1,b3,D12) #gdp_ij(ypc2,ypc1).*a - pop1.*b #gpd as plus and population size as minus 
m21 = logit(expo21)#logit(b1,gdp2,b2,0,b3,0)# logit(b1,gdp2,b2,pop2,b3,D21)#gdp_ij(ypc1,ypc2).*a - pop2.*b
m21_gpc = logit(expo21_gpc)

Mij(P,m) = P.*m
@register Mij(P,m)

M12= Mij(P₁,0)
M21= Mij(P₂,m21)
M21_gpc= Mij(P₂_gpc,m21_gpc)


# --> system looses population
migallow = 1

#differential equations

eqs = [
        D(P₁) ~ r₁*P₁*(1 - P₁/Kp₁) - M12 + M21*migallow, #human population dynamics 
        D(P₂) ~ r₂*P₂*(1 - P₂/Kp₂) + M12 - M21, # standard logistic.
        D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - σ₁*exp(G-G₁)*K₁, #change in capital stock = 
        D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - σ₂*exp(G-G₁)*K₂, # savings-depreciation-damages
        D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8) + α*θ(G-G₀), # global externality.
        D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
        D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
        D(e₂) ~ -z*re₂*e₂, # decarbonization in region 2
        D(P_kp1) ~ r₂*P_kp1*(1 - P_kp1/kp_1),
        D(P_kp) ~ r₂*P_kp*(1 - P_kp/kp_2), ## reference for carrying capacity driver 
    ########
    D(P₁_gpc) ~ r₁*P₁_gpc*(1 - P₁_gpc/Kp₁) + M21_gpc*migallow, #human population dynamics 
    D(P₂_gpc) ~ r₂*P₂_gpc*(1 - P₂_gpc/Kp₂) - M21_gpc, # standard logistic.
    D(K₁_gpc) ~ s₁*Yd₁_gpc - δ₁*K₁_gpc - σ₁*exp(G_gpc-G₁)*K₁_gpc, #change in capital stock = 
    D(K₂_gpc) ~ s₂*Yd₂_gpc - δ₂*K₂_gpc - σ₂*exp(G_gpc-G₁)*K₂_gpc, # savings-depreciation-damages
    D(G_gpc)  ~ er₁_gpc*Y₁_gpc + er₂_gpc*Y₂_gpc - u*(G_gpc-2.8) + α*θ(G_gpc-G₀), # global externality.
    D(z_gpc)  ~ η*θ(G_gpc-Tg)*(1-z_gpc), #initiate decarbonization
    D(e₁_gpc) ~ -z_gpc*re₁*e₁_gpc, # decarbonization in region 1
    D(e₂_gpc) ~ -z_gpc*re₂*e₂_gpc, 
    ### extending by variables without migration 
    D(P1_nomig) ~ r₁*P1_nomig*(1 - P1_nomig/Kp₁),
    D(P2_nomig) ~ r₂*P2_nomig*(1 - P2_nomig/Kp₂),
    D(K1_nomig) ~ s₁*Yd₁_nomig - δ₁*K1_nomig - σ₁*exp(G_nomig-G₁)*K1_nomig,
    D(K2_nomig) ~ s₂*Yd₂_nomig - δ₂*K2_nomig - σ₂*exp(G_nomig-G₁)*K2_nomig,
    D(G_nomig) ~ er₁_nomig*Y₁_nomig + er₂_nomig*Y₂_nomig - u*(G_nomig-2.8) + α*θ(G_nomig-G₀),
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
  
  #per capita gni
  plot(sol[Y₁_nomig/P1_nomig],color = "dodgerblue2", label = "HICs",xticks = 0:200:1000);
  plot!(sol[Y₂_nomig/P2_nomig],color = "orangered", label = "LICs", ylims= (0,40));
  
  # plot!(sol[Y₁_gp/P₁_gp],color = "dodgerblue2", label = false,linestyle=:dashdotdot);
  # plot!(sol[Y₂_gp/P₂_gp],color = "orangered", label = false, ylims= (0,40),linestyle=:dashdotdot);
  
  plot!(sol[Y₁_gpc/P₁_gpc],color = "dodgerblue2", label = false,linestyle=:dashdot);
  plot!(sol[Y₂_gpc/P₂_gpc],color = "orangered", label = false, ylims= (0,40),linestyle=:dashdot);
  
  
  # plot!(sol[Y₁_gpcr/P₁_gpcr],color = "dodgerblue2", label = false,linestyle=:dot);
  # plot!(sol[Y₂_gpcr/P₂_gpcr],color = "orangered", label = false, ylims= (0,40),linestyle=:dot);
  
  plot!(sol[Y₁/P₁],color = "dodgerblue2", label = false,linestyle=:dash);
  plot!(sol[Y₂/P₂],xlabel="Time", ylabel="\$/10³",color = "orangered", label = false, linestyle=:dash,
  title = "Per Capita GNI and Externality", legend_column = -1,legend_position = :topleft,  ylims= (0,40));
  
  plot!(twinx(), sol,idxs=[100*G_nomig],color="palegreen4", label = "CO2",ylims = (0,850),xlabel="",  legend_position = :topright);
  #plot!(twinx(), sol,idxs=[100*G_gp],color="palegreen4", label = false,ylims = (0,850),linestyle=:dashdotdot);
  plot!(twinx(), sol,idxs=[100*G_gpc],color="palegreen4", label = false,ylims = (0,850),xlabel="",linestyle=:dashdot);
  #plot!(twinx(), sol,idxs=[100*G_gpcr], color="palegreen4", label = false,ylims = (0,850),linestyle=:dot);
  
  p1=plot!(twinx(), sol,idxs=[100*G],xlabel="", ylabel="Atm CO2(ppm)", color="palegreen4", label = false,
 legend_column = -1, ylims = (0,850),linestyle=:dash);
  
  #population
  
  plot(sol[P1_nomig],color = "dodgerblue2", label = "HICs");
  plot!(sol[P2_nomig],color = "orangered", label = "LICs");
  plot!(sol[P2_nomig+P1_nomig],xlabel=false, color="grey", label = "World");
  
  #plot!(sol[P₁_gp],color = "dodgerblue2", label = false,linestyle=:dashdotdot);
  #plot!(sol[P₂_gp],color = "orangered", label = false,linestyle=:dashdotdot);
  #plot!(sol[P₁_gp+P₂_gp],color="grey", label = false,linestyle=:dashdotdot);
  
  plot!(sol[P₁_gpc],color = "dodgerblue2", label = false,linestyle=:dashdot);
  plot!(sol[P₂_gpc],color = "orangered", label = false, ylims= (0,40),linestyle=:dashdot);
  plot!(sol[P₂_gpc+P₁_gpc],color="grey", label = false,linestyle=:dashdot);
  
  #plot!(sol[P₁_gpcr],color = "dodgerblue2", label = false,linestyle=:dot);
  #plot!(sol[P₂_gpcr],color = "orangered", label = false, ylims= (0,40),linestyle=:dot);
  #plot!(sol[P₂_gpcr+P₁_gpcr], color="grey", label = false,linestyle=:dot);
  
  plot!(sol,idxs=[P₁],color = "dodgerblue2", label = false ,linestyle=:dash);
  plot!(sol,idxs=[P₂],color = "orangered", label = false,linestyle=:dash);
  p2=plot!(sol[P₁+P₂],xlabel="Time", ylabel="Billion People", label=false, 
  title = "Population",legend_column = -1,legend_position = :topleft, yticks = 0:5:20,
  color = "gray", ylims= (0,15),linestyle=:dash); 
  
  #Distribution 
  plot(x,y,color = "orange", label = false);
  plot!(sol[Y₁_nomig/P1_nomig],sol[Y₂_nomig/P2_nomig], color = "black", label = false);
  #plot!(sol[Y₁_gp/P₁_gp],sol[Y₂_gp/P₂_gp], color = "black", label = false,linestyle=:dashdotdot);
  plot!(sol[Y₁_gpc/P₁_gpc],sol[Y₂_gpc/P₂_gpc], color = "black", label = false,linestyle=:dashdot);
  #plot!(sol[Y₁_gpcr/P₁_gpcr],sol[Y₂_gpcr/P₂_gpcr], color = "black", label = false,linestyle=:dot);
  p3=plot!(sol[Y₁/P₁],sol[Y₂/P₂], xlabel="Per Capita GNI, HICs (\$/10³)",linestyle=:dash,
  ylabel="Per Capita GNI, LICs (\$/10³)", color = "black", label = false, title = "JOS");
  
  
  #operating space
  plot(sol[G_nomig*100],sol[Y₁_nomig+Y₂_nomig], color= "black", label = false);
  # plot!(sol[G_gp*100],sol[Y₁_gp+Y₂_gp], color= "black", label = false,linestyle=:dashdotdot);
  plot!(sol[G_gpc*100],sol[Y₁_gpc+Y₂_gpc], color= "black", label = false,linestyle=:dashdotdot);
  #plot!(sol[G_gpcr*100],sol[Y₁_gpcr+Y₂_gpcr], color= "black", label = false,linestyle=:dashdotdot);
  p4=plot!(sol[G*100],sol[Y₁+Y₂], xlabel="Atm CO2(ppm)",
  ylabel="Total World GNI(\$/10¹²)", color= "black", label = false, title="SOS",linestyle=:dash);
  
  plt = plot(p2, p1, p3, p4, layout=(2,2), size=(750,500))
  display(plt)


#safe

#savefig(plt,"test.png")

#colors: 1(HIC) color = "dodgerblue2", 2 (LIC)color = "orangered", CO2  color="palegreen4"


x = LinRange(1,1,1000)
plot((sol[(Y₁*P₂)/(P₁*Y₂)]),label="a=0.04,b=0.04,r=0")
#plot!((sol[(Y₁_gp*P₂_gp)/(P₁_gp*Y₂_gp)]),label="a=0.04,b=0.03,r=0")
plot!((sol[(Y₁_gpc*P₂_gpc)/(P₁_gpc*Y₂_gpc)]),label="a=0.03,b=0.04,r=0")
#plot!((sol[(Y₁_gpcr*P₂_gpcr)/(P₁_gpcr*Y₂_gpcr)]),label="a=0.03,b=0.03,r=0")
#p5 = plot!((sol[Y₁_nomig/P1_nomig]/sol[Y₂_nomig/P2_nomig]), title = "PC1/PC2")
#plot!(x,label="just",color=:red)
p5 = plot!((sol[(Y₁_nomig*P2_nomig)/ (P1_nomig*Y₂_nomig)]), ylabel = "YPC1/YPC2", xlabel = "Time", title="Per capita Income ratio (ypc1/ypc2)", label = "m=0")

plot(sol[G*100],(sol[(Y₁*P₂)/(P₁*Y₂)]),label="a=0.04,b=0.04,r=0")
#plot!(sol[G_gp*100],(sol[(Y₁_gp*P₂_gp)/(P₁_gp*Y₂_gp)]),label="a=0.04,b=0.03,r=0")
plot!(sol[G_gpc*100],(sol[(Y₁_gpc*P₂_gpc)/(P₁_gpc*Y₂_gpc)]),label="a=0.03,b=0.04,r=0")
#plot!(sol[G_gpcr*100],(sol[(Y₁_gpcr*P₂_gpcr)/(P₁_gpcr*Y₂_gpcr)]),label="a=0.03,b=0.03,r=0")
#p5 = plot!((sol[Y₁_nomig/P1_nomig]/sol[Y₂_nomig/P2_nomig]), title = "PC1/PC2")
#plot!(x,label="just",color=:red)
p6 = plot!(sol[G_nomig*100],(sol[(Y₁_nomig*P2_nomig)/ (P1_nomig*Y₂_nomig)]), ylabel = "YPC1/YPC2", xlabel = "Atm CO2 (ppm)", title="Per capita Income (ypc1/ypc2) \n Carbon Concentration", label = "m=0")


plot(sol[P₁/P₂],label="a=0.04, b=0.04, r=0")
#plot!(sol[P₁_gp/P₂_gp],label="a=0.04, b=0.03, r=0")
plot!(sol[P₁_gpc/P₂_gpc],label="a=0.03,b=0.04,r=0")
#plot!(sol[P₁_gpcr/P₂_gpcr],label="a=0.03,b=0.04,r=0")
p7 = plot!(sol[P1_nomig/P2_nomig], label = "m=0", title = "Popuation ration (P1/P2)", xlabel = "Time", ylabel = "P1/P2")


plot(sol[G*100],label="a=0.04, b=0.04, r=0")
#plot!(sol[G_gp*100],sol[P₁_gp./P₂_gp],label="a=0.04, b=0.03, r=0")
plot!(sol[G_gpc*100],label="a=0.03,b=0.04,r=0")
#plot!(sol[G_gpcr*100],sol[P₁_gpcr./P₂_gpcr],label="a=0.03,b=0.04,r=0")
p8 = plot!(sol[G_nomig]*100, label = "m=0", title = "and Atm CO2", xlabel = "", ylabel = "atm CO2 (ppm)")


plt2 = plot(p5,p6,p7,p8 ,layout = (2,2),size=(900,900))
display(plt2)