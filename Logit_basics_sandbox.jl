
using ModelingToolkit, DifferentialEquations, Plots


@parameters t
D = Differential(t)

r_vary=LinRange(0,0.15,4)

parspop = @parameters r₁=0.038 Kp₁=1.5 r₂=0.042 Kp₂=9.7 mₒ=5 m₁=5 #population/migration dynamics
parsmig = @parameters r = 0 b1= 0.02 b2= 0.02 b3= 0.0  a = 0#migration parameter
parsecn = @parameters δ₁=0.05 δ₂=0.05 α₁=0.5 α₂=0.5 s₁=0.25 s₂=0.21 Cₘ₁=0.7 Cₘ₂=0.7 a₁=2.7 a₂=1.7 #economic dynamics
parsear = @parameters re₁=0.1 re₂=0.1 eb=0.00004 u=0.0025 α=0.1 σ₁=0.03 σ₂=0.03 G₀=20 G₁=5 Tg=50 η=1 #earth system and climate 
pars = [parspop;parsmig;parsecn;parsear] #combine into a single vector.

vars = @variables P_kp_gpcr(t) = 0.24 P_kp_gpc(t) = 0.24  P_kp_gp(t) = 0.24  P₁(t)=0.24 P₂(t)=0.24 K₁(t)=0 K₂(t)=0 G(t)=2.8 e₁(t)=0.0004 e₂(t)=0.0004 z(t)=0 P1_nomig(t)=0.24 P2_nomig(t)=0.24 K1_nomig(t)=0 K2_nomig(t)=0 G_nomig(t)=2.8 e1_nomig(t)=0.0004 e2_nomig(t)=0.0004 z_nomig(t)=0 P₁_gp(t)=0.24 P₂_gp(t)=0.24 K₁_gp(t)=0 K₂_gp(t)=0 G_gp(t)=2.8 e₁_gp(t)=0.0004 e₂_gp(t)=0.0004 z_gp(t)=0 P₁_gpc(t)=0.24 P₂_gpc(t)=0.24 K₁_gpc(t)=0 K₂_gpc(t)=0 G_gpc(t)=2.8 e₁_gpc(t)=0.0004 e₂_gpc(t)=0.0004 z_gpc(t)=0 P₁_gpcr(t)=0.24 P₂_gpcr(t)=0.24 P_kp(t)=0.24 K₁_gpcr(t)=0 K₂_gpcr(t)=0 G_gpcr(t)=2.8 e₁_gpcr(t)=0.0004 e₂_gpcr(t)=0.0004 z_gpcr(t)=0

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
Yd₁_gpcr = Y(P₁_gpcr,Yᵢ₁_gpcr)
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
D_21_g = dia(P₁,P1_nomig)
D_21_gp = dia(P₁_gp,P1_nomig)
D_21_gpc = dia(P₁_gpc,P1_nomig)
D_21_gpcr = dia(P₁_gpcr,P1_nomig)

# Yd₁ = ϕ(Y₁ - P₁*Cₘ₁ - r*(D_21_g*(Y₁/P₁)))
# Yd₂ = ϕ(Y₂ - P₂*Cₘ₂ + r*(D_21_g*(Y₁/P₁)))
# Yd₁_gp = ϕ(Y₁_gp - P₁_gp*Cₘ₁+r*(D_21_gp*(Y₁_gp/P₁_gp)))
# Yd₂_gp = ϕ(Y₂_gp - P₂_gp*Cₘ₂-r*(D_21_gp*(Y₁_gp/P₁_gp)))
# Yd₁_gpc = ϕ(Y₁_gpc - P₁_gpc*Cₘ₁ - r*(D_21_gpc*(Y₁_gpc/P₁_gpc)))
# Yd₂_gpc = ϕ(Y₂_gpc - P₂_gpc*Cₘ₂ + r*(D_21_gpc*(Y₁_gpc/P₁_gpc)))
# Yd₁_gpcr = ϕ(Y₁_gpcr - P₁_gpcr*Cₘ₁ - r*(D_21_gpcr*(Y₁_gpcr/P₁_gpcr)))
# Yd₂_gpcr = ϕ(Y₂_gpcr - P₂_gpcr*Cₘ₂ + r*(D_21_gpcr*(Y₁_gpcr/P₁_gpcr)))
# Yd₁_nomig = ϕ(Y₁_nomig - P1_nomig*Cₘ₁)
# Yd₂_nomig = ϕ(Y₂_nomig - P2_nomig*Cₘ₂)


Yd₁ = ϕ(Y₁ - P₁*Cₘ₁ - r_vary[1]*(D_21_g*(Y₁/P₁)))
Yd₂ = ϕ(Y₂ - P₂*Cₘ₂ + r_vary[1]*(D_21_g*(Y₁/P₁)))
Yd₁_gp = ϕ(Y₁_gp - P₁_gp*Cₘ₁+r_vary[2]*(D_21_gp*(Y₁_gp/P₁_gp)))
Yd₂_gp = ϕ(Y₂_gp - P₂_gp*Cₘ₂-r_vary[2]*(D_21_gp*(Y₁_gp/P₁_gp)))
Yd₁_gpc = ϕ(Y₁_gpc - P₁_gpc*Cₘ₁ - r_vary[3]*(D_21_gpc*(Y₁_gpc/P₁_gpc)))
Yd₂_gpc = ϕ(Y₂_gpc - P₂_gpc*Cₘ₂ + r_vary[3]*(D_21_gpc*(Y₁_gpc/P₁_gpc)))
Yd₁_gpcr = ϕ(Y₁_gpcr - P₁_gpcr*Cₘ₁ - r_vary[4]*(D_21_gpcr*(Y₁_gpcr/P₁_gpcr)))
Yd₂_gpcr = ϕ(Y₂_gpcr - P₂_gpcr*Cₘ₂ + r_vary[4]*(D_21_gpcr*(Y₁_gpcr/P₁_gpcr)))
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
pop2_gp = pop_i(P₂_gp,P₁_gp)
pop2_gpc = pop_i(P₂_gpc,P₁_gpc)
pop2_gpcr = pop_i(P₂_gpcr,P₁_gpcr)

#GDP difference function (share of per capita at home vs per capita destination)
ypc1 = Y₁./P₁
ypc2 = Y₂./P₂
ypc1_gp = Y₁_gp./P₁_gp
ypc2_gp = Y₂_gp./P₂_gp
ypc1_gpc = Y₁_gpc./P₁_gpc
ypc2_gpc = Y₂_gpc./P₂_gpc
ypc1_gpcr = Y₁_gpcr./P₁_gpcr
ypc2_gpcr = Y₂_gpcr./P₂_gpcr

## from more people to poor to move
gdp_ij(ypc1,ypc2) = (ypc1/(ypc2+ypc1))#.*(1/(1+exp.(γ*((Yi./Pi)-Ypc0)))) ##(ypc1/(ypc2+ypc1))
@register gdp_ij(ypc1,ypc2)

gdp1 = gdp_ij(ypc1,ypc2)
gdp2 = gdp_ij(ypc2,ypc1)
gdp2_gp = gdp_ij(ypc2_gp,ypc1_gp)
gdp2_gpc = gdp_ij(ypc2_gpc,ypc1_gpc)
gdp2_gpcr = gdp_ij(ypc2_gpcr,ypc1_gpcr)

##depr Carrying capacity for 
kp(Kp,G) = Kp - G*θ(G-4)#(ϕ(G-4))
@register kp(Kp,G)
kp_2 = kp(Kp₂,G)
kp_2_gp = kp(Kp₂,G_gp)
kp_2_gpc = kp(Kp₂,G_gpc)
kp_2_gpcr = kp(Kp₂,G_gpcr)

##constant
expo21 = expo(4,gdp2,-3,pop2,0.5,θ(P₂-P_kp))
expo21_gp = expo(4,gdp2_gp,-3,pop2_gp,0.5,θ(P₂-P_kp_gp))
expo21_gpc = expo(4,gdp2_gpc,-3,pop2_gpc,0.5,θ(P₂-P_kp_gpcr))
expo21_gpcr = expo(4,gdp2_gpcr,-3,pop2_gpcr,0.5,θ(P₂-P_kp_gpcr))
expo12 = 0 #expo(1,0,-1,0,0.1,D21)

###vary 
# expo21 = expo(4,gdp2,-3,pop2,0.5,θ(P₂-P_kp))
# expo21_gp = expo(4,gdp2_gp,-2,pop2_gp,0.5,θ(P₂-P_kp_gp))
# expo21_gpc = expo(4,gdp2_gpc,-1,pop2_gpc,0.5,θ(P₂-P_kp_gpcr))
# expo21_gpcr = expo(3,gdp2_gpcr,-3,pop2_gpcr,0.5,θ(P₂-P_kp_gpcr))
# expo12 = 0 #expo(1,0,-1,0,0.1,D21)

#pij(Yj,Pj,Yi,Pi) = 1 .-(Yi./(Pi.*((Yi./Pi)+(Yj./Pj)))) 
#pij(Yj,Pj,Yi,Pi) = 1 -(Yj/(Pj*((Yj/Pj)+(Yi/Pi)))) ## falschrum
#@register pij(Yj,Pj,Yi,Pi)

#logit(expo12)#logit(b1,gdp1,b2,pop1,b3,D12) #gdp_ij(ypc2,ypc1).*a - pop1.*b #gpd as plus and population size as minus 
m21 = logit(expo21)#logit(b1,gdp2,b2,0,b3,0)# logit(b1,gdp2,b2,pop2,b3,D21)#gdp_ij(ypc1,ypc2).*a - pop2.*b
m21_gp = logit(expo21_gp)
m21_gpc = logit(expo21_gpc)
m21_gpcr = logit(expo21_gpcr)

Mij(P,m) = P.*m
@register Mij(P,m)

M21= Mij(P₂,m21)
M21_gp= Mij(P₂_gp,m21_gp)
M21_gpc= Mij(P₂_gpc,m21_gpc)
M21_gpcr= Mij(P₂_gpcr,m21_gpcr)

# --> system looses population
migallow = 1

#differential equations

eqs = [
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
  D(P_kp_gp) ~ r₂*P_kp_gp*(1 - P_kp/kp_2_gp),
  D(P_kp_gpc) ~ r₂*P_kp_gpc*(1 - P_kp_gpc/kp_2_gpc),
  D(P_kp_gpcr) ~ r₂*P_kp_gpcr*(1 - P_kp_gpcr/kp_2_gpcr),
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
    ]

@named sys = ODESystem(eqs,t,vars,pars);


prob = ODEProblem(sys, [], (0.0,1000), [])
sol = solve(prob,saveat=1.0)
    



##population on time 
plot(sol[P₁/P₂],label="b1=4,b2=-3,b3=0.5,r=0")
plot!(sol[P₁_gp/P₂_gp],label="b1=4,b2=-3,b3=0.5,r=0.05")
plot!(sol[P₁_gpc/P₂_gpc],label="b1=4,b=-3,b3=0.5,r=0.1")
plot!(sol[P₁_gpcr/P₂_gpcr],label="b1=4,b=-3,b3=0.5,r=0.15")
p5 = plot!(sol[P1_nomig/P2_nomig], label = "m=0", title = "(a) \n \n Population Ratio (P1/P2)", xlabel = "Time", ylabel = "P1/P2")


## ypc distribution on time 
plot((sol[(Y₁*P₂)/(P₁*Y₂)]),label=false)
plot!((sol[(Y₁_gp*P₂_gp)/(P₁_gp*Y₂_gp)]),label=false)
plot!((sol[(Y₁_gpc*P₂_gpc)/(P₁_gpc*Y₂_gpc)]),label=false)
plot!((sol[(Y₁_gpcr*P₂_gpcr)/(P₁_gpcr*Y₂_gpcr)]),label=false)
#p5 = plot!((sol[Y₁_nomig/P1_nomig]/sol[Y₂_nomig/P2_nomig]), title = "PC1/PC2")
#plot!(x,label="just",color=:red)
p6 = plot!((sol[(Y₁_nomig*P2_nomig)/ (P1_nomig*Y₂_nomig)]), ylabel = "YPC1/YPC2", xlabel = "Time", title="(b) \n \n Per capita Income Ratio (ypc1/ypc2)", label = false)

### atm conc on time 
plot(sol[G*100],label=false)
plot!(sol[G_gp*100],label=false)
plot!(sol[G_gpc*100],label=false)
plot!(sol[G_gpcr*100],label=false)
p7 = plot!(sol[G_nomig]*100,label = false, title = "(c) \n \n Carbon Concentration", xlabel = "Time", ylabel = "Atm CO2 (ppm)")

## YPC ratio on G
plot(sol[G*100],(sol[P₁/(P₂)]),label=false)
plot!(sol[G_gp*100],(sol[(P₁_gp/P₂_gp)]),label=false)
plot!(sol[G_gpc*100],(sol[(P₁_gpc/P₂_gpc)]),label=false)
plot!(sol[G_gpcr*100],(sol[(P₁_gpcr/P₂_gpcr)]),label=false)
#p5 = plot!((sol[Y₁_nomig/P1_nomig]/sol[Y₂_nomig/P2_nomig]), title = "PC1/PC2")
#plot!(x,label="just",color=:red)
p8 = plot!(sol[G_nomig*100],(sol[(P1_nomig/P2_nomig)]), ylabel = "YPC1/YPC2", xlabel = "Atm CO2 (ppm)", title="(d) \n \n Population Ratio (P1/P2) \n Carbon Concentration", label = false)



plt2 = plot(p5,p6,p7,p8 ,layout = (2,2),size=(900,900))
display(plt2)

