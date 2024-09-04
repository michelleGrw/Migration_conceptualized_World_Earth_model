
using ModelingToolkit, DifferentialEquations, Plots


a_vary = LinRange(2,5,4)
l = length(a_vary)
rows = 1001

 # Example value for l; replace with your actual value
# Initialize the matrices for b0
P1_results_b0 = zeros(rows, l)
P2_results_b0 = zeros(rows, l)
P_results_b0 = zeros(rows, l)
Y1_results_b0 = zeros(rows, l)
Y2_results_b0 = zeros(rows, l)
YPC1_results_b0 = zeros(rows, l)
YPC2_results_b0 = zeros(rows, l)
G_results_b0 = zeros(rows, l)

# Initialize the matrices for b1
P1_results_b1 = zeros(rows, l)
P2_results_b1 = zeros(rows, l)
P_results_b1 = zeros(rows, l)
Y1_results_b1 = zeros(rows, l)
Y2_results_b1 = zeros(rows, l)
YPC1_results_b1 = zeros(rows, l)
YPC2_results_b1 = zeros(rows, l)
G_results_b1 = zeros(rows, l)

# Initialize the matrices for b2
P1_results_b2 = zeros(rows, l)
P2_results_b2 = zeros(rows, l)
P_results_b2 = zeros(rows, l)
Y1_results_b2 = zeros(rows, l)
Y2_results_b2 = zeros(rows, l)
YPC1_results_b2 = zeros(rows, l)
YPC2_results_b2 = zeros(rows, l)
G_results_b2 = zeros(rows, l)

# Initialize the matrices for b3
P1_results_b3 = zeros(rows, l)
P2_results_b3 = zeros(rows, l)
P_results_b3 = zeros(rows, l)
Y1_results_b3 = zeros(rows, l)
Y2_results_b3 = zeros(rows, l)
YPC1_results_b3 = zeros(rows, l)
YPC2_results_b3 = zeros(rows, l)
G_results_b3 = zeros(rows, l)


for i in 1:l

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

expo21 = expo(a_vary[i],gdp2,0,pop2,1,kp_gp)
expo21_gpc = expo(a_vary[i],gdp2_gpc,-1,pop2_gpc,1,kp_gpc)
expo12 = 0 #expo(1,0,-1,0,0.1,D21)

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
    

P1_results_b0[:,i] = sol[P₁]
P2_results_b0[:,i] = sol[P₂]
P_results_b0[:,i] = sol[P₁+P₂]
Y1_results_b0[:,i] = sol[Y₁]
Y2_results_b0[:,i] = sol[Y₂]
G_results_b0[:,i] = sol[G]
YPC1_results_b0[:,i] = sol[Y₁/P₁]
YPC2_results_b0[:,i] = sol[Y₂/P₂]

P1_results_b1[:,i] = sol[P₁_gpc]
P2_results_b1[:,i] = sol[P₂_gpc]
P_results_b1[:,i] = sol[P₁_gpc+P₂_gpc]
Y1_results_b1[:,i] = sol[Y₁_gpc]
Y2_results_b1[:,i] = sol[Y₂_gpc]
G_results_b1[:,i] = sol[G_gpc]
YPC1_results_b1[:,i] = sol[Y₁_gpc/P₁_gpc]
YPC2_results_b1[:,i] = sol[Y₂_gpc/P₂_gpc]


end



for i in 1:l

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
  
  expo21 = expo(a_vary[i],gdp2,-2,pop2,1,kp_gp)
  expo21_gpc = expo(a_vary[i],gdp2_gpc,-3,pop2_gpc,1,kp_gpc)
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
      
  
  P1_results_b2[:,i] = sol[P₁]
  P2_results_b2[:,i] = sol[P₂]
  P_results_b2[:,i] = sol[P₁+P₂]
  Y1_results_b2[:,i] = sol[Y₁]
  Y2_results_b2[:,i] = sol[Y₂]
  G_results_b2[:,i] = sol[G]
  YPC1_results_b2[:,i] = sol[Y₁/P₁]
  YPC2_results_b2[:,i] = sol[Y₂/P₂]
  
  P1_results_b3[:,i] = sol[P₁_gpc]
  P2_results_b3[:,i] = sol[P₂_gpc]
  P_results_b3[:,i] = sol[P₁_gpc+P₂_gpc]
  Y1_results_b3[:,i] = sol[Y₁_gpc]
  Y2_results_b3[:,i] = sol[Y₂_gpc]
  G_results_b3[:,i] = sol[G_gpc]
  YPC1_results_b3[:,i] = sol[Y₁_gpc/P₁_gpc]
  YPC2_results_b3[:,i] = sol[Y₂_gpc/P₂_gpc]
  
  
  end
  

  
  p1 = plot(ylabel="b2=0 \n Billion People", title = "b1=2",label=false, ylim=(0, 12))
  p2 = plot(title = "b1=3",label=false, ylim=(0, 12))
  p3 = plot(title = "b1=4",label=false, ylim=(0, 12))
  p4 = plot(title = "b1=5",label=false, ylim=(0, 12))
  p5 = plot(ylabel="b2=-1 \n Billion People", label=false, ylim=(0, 12))
  p6 = plot(label=false, ylim=(0, 12))
  p7 = plot(label=false, ylim=(0, 12))
  p8 = plot(label=false, ylim=(0, 12))
  p9 = plot(ylabel="b2=-2 \n Billion People", label=false, ylim=(0, 12))
  p10 = plot(label=false, ylim=(0, 12))
  p11 = plot(label=false, ylim=(0, 12))
  p12 = plot(label=false, ylim=(0, 12))
  p13 = plot(ylabel="b2=-3 \n Billion People", label=false, ylim=(0, 12))
  p14 = plot(label=false, ylim=(0, 12))
  p15 = plot(label=false, ylim=(0, 12))
  p16 = plot(label=false, ylim=(0, 12))
 



#### for b0 
#a1
    plot!(p1, P1_results_b0[:, 1], color="dodgerblue2",label=false)
    plot!(p1, P2_results_b0[:, 1], color="orangered", label=false)
    plot!(p1, P_results_b0[:, 1], color="gray", label=false)
#a2
    plot!(p2, P1_results_b0[:, 2], color="dodgerblue2",label=false)
    plot!(p2, P2_results_b0[:, 2], color="orangered", label=false)
    plot!(p2, P_results_b0[:, 2], color="gray", label=false)
#a3
plot!(p3, P1_results_b0[:, 3], color="dodgerblue2",label=false)
plot!(p3, P2_results_b0[:, 3], color="orangered", label=false)
plot!(p3, P_results_b0[:, 3], color="gray", label=false)
#a4
plot!(p4, P1_results_b0[:, 4], color="dodgerblue2",label=false)
plot!(p4, P2_results_b0[:, 4], color="orangered", label=false)
plot!(p4, P_results_b0[:, 4], color="gray",label=false)

###b1
#a1
plot!(p5, P1_results_b1[:, 1], color="dodgerblue2",label=false)
plot!(p5, P2_results_b1[:, 1], color="orangered", label=false)
plot!(p5, P_results_b1[:, 1], color="gray", label=false)
#a2
plot!(p6, P1_results_b1[:, 2], color="dodgerblue2",label=false)
plot!(p6, P2_results_b1[:, 2], color="orangered", label=false)
plot!(p6, P_results_b1[:, 2], color="gray", label=false)
#a3
plot!(p7, P1_results_b1[:, 3], color="dodgerblue2",label=false)
plot!(p7, P2_results_b1[:, 3], color="orangered", label=false)
plot!(p7, P_results_b1[:, 3], color="gray", label=false)
#a4
plot!(p8, P1_results_b1[:, 4], color="dodgerblue2",label=false)
plot!(p8, P2_results_b1[:, 4], color="orangered", label=false)
plot!(p8, P_results_b1[:, 4], color="gray", label=false)

#### b2
#a1
plot!(p9, P1_results_b2[:, 1], color="dodgerblue2",label=false)
plot!(p9, P2_results_b2[:, 1], color="orangered",label=false)
plot!(p9, P_results_b2[:, 1], color="gray", label=false)
#a2
plot!(p10, P1_results_b2[:, 2], color="dodgerblue2",label=false)
plot!(p10, P2_results_b2[:, 2], color="orangered", label=false)
plot!(p10, P_results_b2[:, 2], color="gray", label=false)
#a3
plot!(p11, P1_results_b2[:, 3], color="dodgerblue2",label=false)
plot!(p11, P2_results_b2[:, 3], color="orangered", label=false)
plot!(p11, P_results_b2[:, 3], color="gray", label=false)
#a4
plot!(p12, P1_results_b2[:, 4], color="dodgerblue2",label=false)
plot!(p12, P2_results_b2[:, 4], color="orangered", label=false)
plot!(p12, P_results_b2[:, 4], color="gray", label=false)
###b3
#a1
plot!(p13, P1_results_b3[:, 1], color="dodgerblue2",label=false)
plot!(p13, P2_results_b3[:, 1], color="orangered", label=false)
plot!(p13, P_results_b3[:, 1], color="gray", label=false)
#a2
plot!(p14, P1_results_b3[:, 2], color="dodgerblue2",label=false)
plot!(p14, P2_results_b3[:, 2], color="orangered", label=false)
plot!(p14, P_results_b3[:, 2], color="gray",label=false)
#a3
plot!(p15, P1_results_b3[:, 3], color="dodgerblue2",label=false)
plot!(p15, P2_results_b3[:, 3], color="orangered", label=false)
plot!(p15, P_results_b3[:, 3], color="gray",label=false)
#a4
plot!(p16, P1_results_b3[:, 4], color="dodgerblue2",label=false)
plot!(p16, P2_results_b3[:, 4], color="orangered",label=false)
plot!(p16, P_results_b3[:, 4], color="gray", label=false)

plt = plot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,layout=(4,4),size=(900,900))
display(plt)


############### relations 



# # # ###GNI per capita 
# p2=plot(xlabel="Time", ylabel="\$/10³", title = "Per capita GNI",label = false)#ylabel="Atmospheric CO2(ppm)",yticks = 0:10:75, ylims = (0,85));
# for i in 1:l
# plot!(p2,YPC1_results_b0[:,i],color="dodgerblue2",linestyle=:dot, label = false);
# plot!(p2,YPC2_results_b0[:,i],color = "orangered", linestyle=:dash,label = false)
# end

# ##atm co2
# p3=plot(xlabel="Time", ylabel="ppm", title = "Atm Carbon Concentration",label =false)#ylabel="Atmospheric CO2(ppm)",yticks = 0:10:75, ylims = (0,85));
# for i in 1:l
#     a = (i-1)*0.01
# plot!(p3,G_results_b0[:,i],label="$a")
# end

# # # #Distribution
# p4=plot( xlabel="per capita GNI, HICs (\$/10³)", ylabel="per capita GNI, LICs (\$/10³)", title = "GNI distribution")
# for i in 1:l
#     a = (i-1)*0.01
#     plot!(p4,YPC1_results_b0[:,i],YPC2_results_b0[:,i], label = "$a",xlims = (0,40));
#     end
# plot!(p4,x,y,color = "red", label ="just")

# # # #operating space
# p5 = plot(xlabel="Atmospheric CO2(ppm)",ylabel="Total World GNI(\$/10¹²)", title="Operating Space")
# for i in 1:l
#     a = (i-1)*0.01
# plot!(p5,G_results_b0[:,i]*100,Y1_results_b0[:,i].+Y2_results_b0[:,i],label = "$a")
# end


# plt = plot(p1,p2,p3,p4,p5,layout=(3,2),size=(900,700))
