# ---------- #
# JUnit tests expected values
# ---------- #

library(ape)
library(diversitree)
library(hisse)

# ----- START: BiSSE validation ----- #
# JUnit: HSDSEPBiSSETest
trstr <- "(((sp15:10.27880339,(sp57:0.4327353378,sp58:0.4327353378):9.846068053):21.30935137,((((sp49:1.322566942,sp50:1.322566942):6.531246386,(((((sp42:1.618558172,sp43:1.618558172):1.249323508,sp37:2.86788168):0.4105311845,sp36:3.278412865):1.110829025,sp28:4.38924189):2.453996398,((sp53:0.6765630317,sp54:0.6765630317):5.834067793,sp21:6.510630824):0.3326074635):1.01057504):6.546385565,sp12:14.40019889):3.891878236,((((sp18:8.595427361,((sp19:6.988162304,((sp39:1.941330272,(sp59:0.4256083779,sp60:0.4256083779):1.515721894):1.374985348,sp35:3.31631562):3.671846684):1.028692949,(sp24:5.527011086,(sp25:5.478875203,(sp40:1.898502308,sp41:1.898502308):3.580372894):0.04813588287):2.489844168):0.5785721075):0.8605508177,((sp47:1.324188282,sp48:1.324188282):1.210143714,sp38:2.534331996):6.921646183):1.848794077,(sp22:6.144323416,sp23:6.144323416):5.160448839):4.752352041,sp10:16.0571243):2.234952832):13.29607763):8.9940146,(sp6:33.80408947,(((sp29:4.271294196,sp30:4.271294196):3.963360008,(sp46:1.515605972,(sp51:0.6842469553,sp52:0.6842469553):0.8313590168):6.719048232):21.69107479,((((sp44:1.517683119,sp45:1.517683119):13.83340518,((sp33:3.451233406,sp34:3.451233406):7.318030201,sp14:10.76926361):4.581824694):2.3268441,((sp31:3.988873926,sp32:3.988873926):13.39833,(sp26:5.46221229,sp27:5.46221229):11.92499164):0.2907284735):12.10203097,((sp16:9.676541191,sp17:9.676541191):11.55054389,(sp11:16.00734921,(sp55:0.6152478573,sp56:0.6152478573):15.39210136):5.219735869):8.552878292):0.1457656227):3.878360468):6.778079891);"
tr = read.tree(file="", text=trstr)
states <- c(1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2) - 1
names(states) <- tr$tip.label

sampling.f <- c(1,1)

bisselk <- make.bisse(tree=tr, states=states, sampling.f=sampling.f, strict=FALSE)

bisse.pars <- c(0.15, 0.3, 0.1, 0.1, 0.05, 0.05)
names(bisse.pars) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")

res <- bisselk(pars=bisse.pars, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE) # prior on root is the weighted average of D0 and D1, i.e., ROOT.OBS = D = D0 * (D0/(D0+D1)) + D1 * (D1/(D0+D1))

res[1] # -198.2515
# ----- END: BiSSE validation ----- #

# ----- START: MuSSE validation ----- #
# JUnit: HSDSEPMuSSETest
# pars for inference
pars <- c(.1,  .15,  .2, .1, # lambda 1, 2, 3, 4
.03, .045, .06, 0.03, # mu 1, 2, 3, 4
.05, .05, .00, # q12, q13, q14
.05, .00, .05, # q21, q23, q24
.05, .00, .05, # q31, q32, q34
.00, .05, .05)

set.seed(2)
phy <- tree.musse(pars, 30, x0=1)
sim.dat <- phy$tip.state

lik <- make.musse(phy, states, 4, strict=FALSE, sampling.f=c(1,1,1,1))
diversitree.free = lik(pars, root.p=NULL, root=ROOT.FLAT, condition.surv=FALSE, intermediates=TRUE)
diversitree.free[1] # -122.8801
# ----- END: MuSSE validation ----- #

# ----- START: MuHiSSE validation ----- #
# JUnit: HSDSEPMuHiSSETest (but same likelihood as above)
# using simulation from above
states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
states <- states[phy$tip.label,]

# recoding states
states.trans <- states
for (i in 1:Ntip(phy)){
    if (states[i,1]==1){
        states.trans[i,1] = 0
        states.trans[i,2] = 0
    }
    if (states[i,1]==2){
        states.trans[i,1] = 0
        states.trans[i,2] = 1
    }
    if (states[i,1]==3){
        states.trans[i,1] = 1
        states.trans[i,2] = 0
    }
    if (states[i,1]==4){
        states.trans[i,1] = 1
        states.trans[i,2] = 1
    }
}

pars.hisse <- rep(c(pars[1]+pars[5],pars[2]+pars[6],pars[3]+pars[7],pars[4]+pars[8],
                pars[5]/pars[1],pars[6]/pars[2],pars[7]/pars[3],pars[8]/pars[4],
                0.05,0.05,0, # ignoring diagonal (so 3 instead of 4)
                0.05,0,0.05, # same for all rows below
                0.05,0,0.05, # same
                0,0.05,0.05, # same
                1, # 00A -> 00B
                0, 0, 0, 0, 0, 0, # 00A -> 00C ... 00H are all 0
                1, # 01A -> 01B
                0, 0, 0, 0, 0, 0,
                1, # 10A -> 10B
                0, 0, 0, 0, 0, 0,
                1, # 11A -> 11B
                0, 0, 0, 0, 0, 0),
                2) # do both top and bottom halves

model.vec <- rep(0,384)
model.vec[1:96] = pars.hisse
phy$node.label = NULL

cache <- hisse:::ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=TRUE,
                                         nb.tip=Ntip(phy), nb.node=Nnode(phy),
                                         bad.likelihood=exp(-500), ode.eps=0)

gen <- hisse:::FindGenerations(phy)
dat.tab <- hisse:::OrganizeData(states.trans, phy, f=c(1,1,1,1), hidden.states=TRUE)
hisse.constrained <- hisse:::DownPassMuHisse(dat.tab, gen=gen, cache=cache, root.p=NULL, root.type="equal", condition.on.survival=FALSE)
round(hisse.constrained,4) # -122.8801
# ----- END: MuHiSSE validation ----- #
           
# ----- START: HiSSE run 1 ----- #
# JUnit: HSDSEPHiSSETest1
set.seed(4)

# Essentially we are setting up a model that models the evolution of two binary characters
# Thus, we are assuming the following state combinations 1=00, 2=10, 3=01, 4=11:
pars <- c(0.1,0.1,0.1,0.2, # lambdas
          rep(0.03,4), # mus
          0.01,0.01,0, # qs
          0.01,0,0.01,
          0.01,0,0.01,
          0,0.01,0.01)

phy <- tree.musse(pars, max.taxa=50, x0=1, include.extinct=FALSE)
write.tree(phy) # for java

sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)
sim.dat[sim.dat[,2]==3,2] = 1 # 3 -> 1
sim.dat[sim.dat[,2]==4,2] = 2 # 4 -> 2
paste0(sim.dat[,1], collapse=",") # for java
paste0(paste(sim.dat[,1], sim.dat[,2], sep="="), collapse=",") # for java

sim.dat[,2] = sim.dat[,2] - 1 # 0|1 instead of 1|2

turnover.anc <- c(1,2,3,4)
eps.anc <- c(1,2,3,4)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates.nodual <- ParDrop(trans.rates,c(3,5,8,10))

set.seed(123)

pp <- hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual, output.type="raw", root.type="equal", condition.on.survival=FALSE, root.p=NULL) # -172.0229

## Fit
##       -lnL      AIC     AICc ntax
##  -172.0229 376.0458 392.5307   50

## Probability of extinction: 

## Diversification Rates
##                          lambda           mu
## rate0A             6.440928e-10 1.417061e-09
## alpha0A            1.000000e+00 1.000000e+00
## beta0A             1.000000e+00 1.000000e+00
## timeslice.factor0A 1.000000e+00 1.000000e+00

##                       lambda         mu
## rate1A             0.1128975 0.04597049
## alpha1A            1.0000000 1.00000000
## beta1A             1.0000000 1.00000000
## timeslice.factor1A 1.0000000 1.00000000

##                       lambda         mu
## rate0B             0.1308483 0.03687651
## alpha0B            1.0000000 1.00000000
## beta0B             1.0000000 1.00000000
## timeslice.factor0B 1.0000000 1.00000000

##                      lambda           mu
## rate1B             85.03492 1.266694e-06
## alpha1B             1.00000 1.000000e+00
## beta1B              1.00000 1.000000e+00
## timeslice.factor1B  1.00000 1.000000e+00

## Transition Rates
##              (0A)        (1A)         (0B)         (1B)
## (0A)           NA 2.11773e-09 2.061154e-09 0.000000e+00
## (1A) 6.432635e-03          NA 0.000000e+00 2.479683e-09
## (0B) 2.355085e-09 0.00000e+00           NA 1.565876e-02
## (1B) 0.000000e+00 1.00000e+02 1.000000e+02           NA

# ----- START: CID2 ----- #
# JUnit: HSDSEPCID2Test
turnover.anc <- c(1,1,2,2)
eps.anc <- c(1,1,2,2)

set.seed(123)

pp <- hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual, output.type="raw", root.type="equal", condition.on.survival=FALSE, root.p=NULL) # -lnL = -178.9531

## Fit
##       -lnL      AIC     AICc ntax
##  -178.9531 381.9063 390.3387   50

## Probability of extinction: 

## Diversification Rates
##                       lambda         mu
## rate0A             0.2108637 0.04285281
## alpha0A            1.0000000 1.00000000
## beta0A             1.0000000 1.00000000
## timeslice.factor0A 1.0000000 1.00000000

##                       lambda         mu
## rate1A             0.2108637 0.04285281
## alpha1A            1.0000000 1.00000000
## beta1A             1.0000000 1.00000000
## timeslice.factor1A 1.0000000 1.00000000

##                       lambda           mu
## rate0B             0.0688838 1.419801e-10
## alpha0B            1.0000000 1.000000e+00
## beta0B             1.0000000 1.000000e+00
## timeslice.factor0B 1.0000000 1.000000e+00

##                       lambda           mu
## rate1B             0.0688838 1.419801e-10
## alpha1B            1.0000000 1.000000e+00
## beta1B             1.0000000 1.000000e+00
## timeslice.factor1B 1.0000000 1.000000e+00

## Transition Rates
##              (0A)         (1A)        (0B)         (1B)
## (0A)           NA 3.054604e-02 0.129056368 0.000000e+00
## (1A) 2.061154e-09           NA 0.000000000 2.001798e-01
## (0B) 5.421624e-02 0.000000e+00          NA 2.061154e-09
## (1B) 0.000000e+00 2.061154e-09 0.008978002           NA

 # ----- END: CID2 ----- #

# ----- START: HiSSE run 2 ----- #
# JUnit: HSDSEPHiSSETest2
# this is the same as run 2, but with hidden trait B having state 1 only (see Fig. 1 in HiSSE paper)
pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
.03, .045, .06, # mu 1, 2, 3
.05, .00, # q12, q13
.05, .05, # q21, q23
.00, .05 # q31, q32
)

set.seed(2)
phy <- tree.musse(pars, 30, x0=1)
write.tree(phy) # for java

sim.dat <- phy$tip.state
sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)
sim.dat[sim.dat[,2]==3,2] = 2 # state 3 -> becomes state 2
paste0(sim.dat[,1], collapse=",") # for java
paste0(paste(sim.dat[,1], sim.dat[,2], sep="="), collapse=",") # for java
sim.dat[,2] = sim.dat[,2] - 1 # 0|1 instead of 1|2

turnover.anc <- c(1,2,0,3)
eps.anc <- c(1,2,0,3)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates.nodual.no0B <- ParDrop(trans.rates,c(2,3,5,7,8,9,10,12))
trans.rates.nodual.no0B

set.seed(123)

pp <- hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.no0B, output.type="raw", root.type="equal", condition.on.survival=FALSE, root.p=NULL) # -lnL = -101.1803

## Fit
##       -lnL      AIC     AICc ntax
##  -101.1803 222.3606 233.9395   30

## Probability of extinction: 

## Diversification Rates
##                        lambda        mu
## rate0A             0.08740614 0.0173183
## alpha0A            1.00000000 1.0000000
## beta0A             1.00000000 1.0000000
## timeslice.factor0A 1.00000000 1.0000000

##                       lambda         mu
## rate1A             0.1505434 0.04541325
## alpha1A            1.0000000 1.00000000
## beta1A             1.0000000 1.00000000
## timeslice.factor1A 1.0000000 1.00000000

##                    lambda mu
## rate0B                  0  0
## alpha0B                 1  1
## beta0B                  1  1
## timeslice.factor0B      1  1

##                      lambda           mu
## rate1B             3.541145 7.427959e-09
## alpha1B            1.000000 1.000000e+00
## beta1B             1.000000 1.000000e+00
## timeslice.factor1B 1.000000 1.000000e+00

## Transition Rates
##            (0A)       (1A) (0B)         (1B)
## (0A)         NA 0.03693569    0 0.000000e+00
## (1A) 0.09483739         NA    0 2.603804e-09
## (0B) 0.00000000 0.00000000   NA 0.000000e+00
## (1B) 0.00000000 4.64567414    0           NA

# ----- START: ClaSSE ----- #
# JUnit: JSDSEPClaSSETest
trstr <- "(((sp39:0.518912972,sp40:0.518912972):19.54206195,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738);"
tr = read.tree(file="", text=trstr)
states <- c(1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 4, 2, 2, 2, 4, 1, 1, 2, 2, 2, 2, 4, 3, 1, 1, 3, 4, 2, 1)
names(states) <- tr$tip.label

sampling.f <- c(1,1,1,1)

classelk = make.classe(tree=tr, states=states, k=4, sampling.f=sampling.f, strict=FALSE)

pars.names = argnames(classelk)
classe.pars = rep(0, times=length(pars.names))
names(classe.pars) = pars.names

# mus
classe.pars[grepl(pattern="mu", x=pars.names)] = 0.1

# qs
classe.pars[pars.names == "q21"] = 0.01
classe.pars[pars.names == "q31"] = 0.01
classe.pars[pars.names == "q24"] = 0.01
classe.pars[pars.names == "q34"] = 0.01
classe.pars[pars.names == "q42"] = 0.01
classe.pars[pars.names == "q43"] = 0.01

# null range
classe.pars[pars.names=="lambda111"] = 0.2

# narrow sympatry
classe.pars[pars.names=="lambda222"] = 0.2
classe.pars[pars.names=="lambda333"] = 0.2

# no jump dispersal
classe.pars[pars.names=="lambda223"] = 0.0
classe.pars[pars.names=="lambda323"] = 0.0

# subsympatry
classe.pars[pars.names=="lambda424"] = 1/6 * 0.2
classe.pars[pars.names=="lambda434"] = 1/6 * 0.2
classe.pars[pars.names=="lambda442"] = 1/6 * 0.2
classe.pars[pars.names=="lambda443"] = 1/6 * 0.2

# allopatry (vicariance)
classe.pars[pars.names=="lambda423"] = 1/6 * 0.2
classe.pars[pars.names=="lambda432"] = 1/6 * 0.2

res <- classelk(pars=classe.pars, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
res[1] # -129.9762
# ----- END: ClaSSE ----- #

# ----- START: ClaHiSSE ----- #
pars.names <- c("lambda111", "lambda112", "lambda113", "lambda114", "lambda115", "lambda116",
"lambda122", "lambda123", "lambda124", "lambda125", "lambda126", "lambda133",
"lambda134", "lambda135", "lambda136", "lambda144", "lambda145", "lambda146",
"lambda155", "lambda156", "lambda166", "lambda211", "lambda212", "lambda213",
"lambda214", "lambda215", "lambda216", "lambda222", "lambda223", "lambda224",
"lambda225", "lambda226", "lambda233", "lambda234", "lambda235", "lambda236",
"lambda244", "lambda245", "lambda246", "lambda255", "lambda256", "lambda266",
"lambda311", "lambda312", "lambda313", "lambda314", "lambda315", "lambda316",
"lambda322", "lambda323", "lambda324", "lambda325", "lambda326", "lambda333",
"lambda334", "lambda335", "lambda336", "lambda344", "lambda345", "lambda346",
"lambda355", "lambda356", "lambda366", "lambda411", "lambda412", "lambda413",
"lambda414", "lambda415", "lambda416", "lambda422", "lambda423", "lambda424",
"lambda425", "lambda426", "lambda433", "lambda434", "lambda435", "lambda436",
"lambda444", "lambda445", "lambda446", "lambda455", "lambda456", "lambda466",
"lambda511", "lambda512", "lambda513", "lambda514", "lambda515", "lambda516",
"lambda522", "lambda523", "lambda524", "lambda525", "lambda526", "lambda533",
"lambda534", "lambda535", "lambda536", "lambda544", "lambda545", "lambda546",
"lambda555", "lambda556", "lambda566", "lambda611", "lambda612", "lambda613",
"lambda614", "lambda615", "lambda616", "lambda622", "lambda623", "lambda624",
"lambda625", "lambda626", "lambda633", "lambda634", "lambda635", "lambda636",
"lambda644", "lambda645", "lambda646", "lambda655", "lambda656", "lambda666",
"mu1", "mu2", "mu3", "mu4", "mu5", "mu6",     
"q12", "q13", "q14", "q15", "q16", "q21",     
"q23", "q24", "q25", "q26", "q31", "q32",
"q34", "q35", "q36", "q41", "q42", "q43",
"q45", "q46", "q51", "q52", "q53", "q54",
"q56", "q61", "q62", "q63", "q64", "q65")

clahisse.pars = rep(0, times=length(pars.names))
names(clahisse.pars) = pars.names

obs.birth.rate <- 0.1
## hidden.birth.rate <- 0.1
hidden.birth.rate <- 0.2

# mus
clahisse.pars[grepl(pattern="mu", x=pars.names)] = 0.01

# lambdas
# observed sympatric
clahisse.pars[pars.names=="lambda111"] = obs.birth.rate
clahisse.pars[pars.names=="lambda222"] = obs.birth.rate
clahisse.pars[pars.names=="lambda333"] = obs.birth.rate

# hidden sympatric
clahisse.pars[pars.names=="lambda444"] = hidden.birth.rate/2
clahisse.pars[pars.names=="lambda555"] = hidden.birth.rate/2
clahisse.pars[pars.names=="lambda666"] = hidden.birth.rate/2

# observed subsympatric
clahisse.pars[pars.names=="lambda313"] = obs.birth.rate
clahisse.pars[pars.names=="lambda323"] = obs.birth.rate

# hidden subsympatric
clahisse.pars[pars.names=="lambda646"] = hidden.birth.rate/2
clahisse.pars[pars.names=="lambda656"] = hidden.birth.rate/2

# observed vicariant (allopatric)
clahisse.pars[pars.names=="lambda312"] = obs.birth.rate

# hidden vicariant (allopatric)
clahisse.pars[pars.names=="lambda645"] = hidden.birth.rate

# qs
clahisse.pars[pars.names == "q12"] = 0.01
clahisse.pars[pars.names == "q13"] = 0.01
clahisse.pars[pars.names == "q14"] = 0.01
clahisse.pars[pars.names == "q15"] = 0.0
clahisse.pars[pars.names == "q16"] = 0.0

clahisse.pars[pars.names == "q21"] = 0.01
clahisse.pars[pars.names == "q23"] = 0.01
clahisse.pars[pars.names == "q24"] = 0.0
clahisse.pars[pars.names == "q25"] = 0.01
clahisse.pars[pars.names == "q26"] = 0.0

clahisse.pars[pars.names == "q31"] = 0.01
clahisse.pars[pars.names == "q32"] = 0.01
clahisse.pars[pars.names == "q34"] = 0.0
clahisse.pars[pars.names == "q35"] = 0.0
clahisse.pars[pars.names == "q36"] = 0.01

clahisse.pars[pars.names == "q41"] = 0.01
clahisse.pars[pars.names == "q42"] = 0.0
clahisse.pars[pars.names == "q43"] = 0.0
clahisse.pars[pars.names == "q45"] = 0.01
clahisse.pars[pars.names == "q46"] = 0.01

clahisse.pars[pars.names == "q51"] = 0.0
clahisse.pars[pars.names == "q52"] = 0.01
clahisse.pars[pars.names == "q53"] = 0.0
clahisse.pars[pars.names == "q54"] = 0.01
clahisse.pars[pars.names == "q56"] = 0.01

clahisse.pars[pars.names == "q61"] = 0.0
clahisse.pars[pars.names == "q62"] = 0.0
clahisse.pars[pars.names == "q63"] = 0.01
clahisse.pars[pars.names == "q64"] = 0.01
clahisse.pars[pars.names == "q65"] = 0.01

set.seed(123)
phy <- tree.classe(clahisse.pars, max.taxa=300, x0=1, include.extinct=FALSE)
sim.dat <- phy$tip.state
names(sim.dat) = names(phy$tip.state)
sim.dat[sim.dat==4] <- 1 # hiding states 4 into 1
sim.dat[sim.dat==5] <- 2 # hiding states 5 into 2
sim.dat[sim.dat==6] <- 3 # hiding states 6 into 3 (widespread into widespread)

## sim.dat[sim.dat==3] <- 0 # recoding it for GeoSSE
## sampling.f <- c(1, 1, 1)
## geosselk <- make.geosse(tree=phy, states=sim.dat, sampling.f=sampling.f, strict=FALSE)

sampling.f <- c(1, 1, 1)

classelk = make.classe(tree=phy, states=sim.dat, k=3, sampling.f=sampling.f, strict=FALSE)
argnames(classelk)
classelk.null = constrain(classelk, lambda111~lambda312, lambda222~lambda312)

classe.pars = c(rep(0.1, 18), rep(0.01, 3), rep(0.01, 6))
fit <- find.mle(classelk, x.init=classe.pars, condition.surv=FALSE)

x <- runif(1, 0.05, 0.2)
y <- runif(1, 0.005, 0.015)
z <- runif(1, 0.005, 0.015)
w <- runif(1, 0.05, 0.2)
classe.pars.null = c(rep(x, 11), w, rep(x, 4), rep(y, 3), rep(z, 6))
fit.null <- find.mle(classelk.null, x.init=classe.pars.null, condition.surv=FALSE)

## res <- clahisselk(pars=clahisse.pars, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
## res[1] # -389.4712
