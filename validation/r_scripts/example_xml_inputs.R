# ---------- #
# Trees and trait data for example .xmls
# ---------- #

library(ape)
library(diversitree)
library(hisse)

# ----- BiSSE_fixed_tree_HSDSEP ----- #
# and
# ----- BiSSE_fixed_tree_SDSEP ----- #
pars <- c(.15, .3, .1, .1, .1, .1) # lambdas, mus, qs

set.seed(12345)
phy <- tree.bisse(pars, max.taxa=60, include.extinct=FALSE, x0=NA)

write.tree(phy)

## "(((sp93:0.2356622822,(sp94:0.2154739279,sp95:0.2154739279):0.02018835424):12.9847595,((((sp57:2.815062535,sp43:2.815062535):1.780249256,(sp36:4.449415911,((((sp99:0.06034332881,sp100:0.06034332881):0.4473268161,sp88:0.5076701449):1.726055497,sp56:2.233725641):1.371042764,((sp72:1.31183414,sp73:1.31183414):1.545878271,sp42:2.857712411):0.7470559947):0.8446475052):0.1458958804):3.207493326,(((sp27:5.033762494,(((((sp71:1.391331894,((sp86:0.5759627483,sp87:0.5759627483):0.7866153889,(sp89:0.4023765821,sp90:0.4023765821):0.9602015551):0.02875375666):1.595833121,(sp77:0.9006605643,sp78:0.9006605643):2.086504451):0.1392348633,sp41:3.126399879):0.8899363982,(((sp48:2.589801526,sp49:2.589801526):0.1590377796,sp46:2.748839306):0.6512535555,(sp60:2.03018895,((sp79:0.8625209504,(sp101:0.007232857647,sp102:0.007232857647):0.8552880928):0.6584553534,sp69:1.520976304):0.5092126459):1.369903912):0.6162434155):0.07444131832,sp34:4.090777595):0.9429848988):0.6018834209,sp25:5.635645915):0.003042611101,((sp59:2.160114586,(sp70:1.408818523,sp81:1.408818523):0.7512960633):2.006513375,(sp50:2.459826639,sp51:2.459826639):1.706801322):1.472060565):2.164116591):1.673667032,(((sp91:0.321558177,sp92:0.321558177):0.6062021719,sp76:0.9277603489):2.437759666,(sp47:2.651955001,(sp54:2.307131853,sp55:2.307131853):0.3448231482):0.7135650134):6.110952134):3.74394963):2.502755808,(sp15:15.70472892,(((sp37:3.930561202,(sp61:1.91883471,sp62:1.91883471):2.011726492):1.624114024,sp24:5.554675226):7.421080767,(((((sp84:0.6179056665,sp85:0.6179056665):0.2681327966,(sp82:0.8310677522,sp83:0.8310677522):0.05497071083):0.6598028498,sp68:1.545841313):3.050734938,sp29:4.596576251):2.433125762,(((sp96:0.1813892509,(sp97:0.1346470644,sp98:0.1346470644):0.04674218648):1.607626581,sp63:1.789015832):1.762305454,(sp39:3.456648649,sp40:3.456648649):0.0946726372):3.478380726):5.946053981):2.728972922):0.01844867133):0.0;"

cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

## <taxon id="sp15" spec="Taxon"/>
## <taxon id="sp24" spec="Taxon"/>
## <taxon id="sp25" spec="Taxon"/>
## <taxon id="sp27" spec="Taxon"/>
## <taxon id="sp29" spec="Taxon"/>
## <taxon id="sp34" spec="Taxon"/>
## <taxon id="sp36" spec="Taxon"/>
## <taxon id="sp37" spec="Taxon"/>
## <taxon id="sp39" spec="Taxon"/>
## <taxon id="sp40" spec="Taxon"/>
## <taxon id="sp41" spec="Taxon"/>
## <taxon id="sp42" spec="Taxon"/>
## <taxon id="sp43" spec="Taxon"/>
## <taxon id="sp46" spec="Taxon"/>
## <taxon id="sp47" spec="Taxon"/>
## <taxon id="sp48" spec="Taxon"/>
## <taxon id="sp49" spec="Taxon"/>
## <taxon id="sp50" spec="Taxon"/>
## <taxon id="sp51" spec="Taxon"/>
## <taxon id="sp54" spec="Taxon"/>
## <taxon id="sp55" spec="Taxon"/>
## <taxon id="sp56" spec="Taxon"/>
## <taxon id="sp57" spec="Taxon"/>
## <taxon id="sp59" spec="Taxon"/>
## <taxon id="sp60" spec="Taxon"/>
## <taxon id="sp61" spec="Taxon"/>
## <taxon id="sp62" spec="Taxon"/>
## <taxon id="sp63" spec="Taxon"/>
## <taxon id="sp68" spec="Taxon"/>
## <taxon id="sp69" spec="Taxon"/>
## <taxon id="sp70" spec="Taxon"/>
## <taxon id="sp71" spec="Taxon"/>
## <taxon id="sp72" spec="Taxon"/>
## <taxon id="sp73" spec="Taxon"/>
## <taxon id="sp76" spec="Taxon"/>
## <taxon id="sp77" spec="Taxon"/>
## <taxon id="sp78" spec="Taxon"/>
## <taxon id="sp79" spec="Taxon"/>
## <taxon id="sp81" spec="Taxon"/>
## <taxon id="sp82" spec="Taxon"/>
## <taxon id="sp83" spec="Taxon"/>
## <taxon id="sp84" spec="Taxon"/>
## <taxon id="sp85" spec="Taxon"/>
## <taxon id="sp86" spec="Taxon"/>
## <taxon id="sp87" spec="Taxon"/>
## <taxon id="sp88" spec="Taxon"/>
## <taxon id="sp89" spec="Taxon"/>
## <taxon id="sp90" spec="Taxon"/>
## <taxon id="sp91" spec="Taxon"/>
## <taxon id="sp92" spec="Taxon"/>
## <taxon id="sp93" spec="Taxon"/>
## <taxon id="sp94" spec="Taxon"/>
## <taxon id="sp95" spec="Taxon"/>
## <taxon id="sp96" spec="Taxon"/>
## <taxon id="sp97" spec="Taxon"/>
## <taxon id="sp98" spec="Taxon"/>
## <taxon id="sp99" spec="Taxon"/>
## <taxon id="sp100" spec="Taxon"/>
## <taxon id="sp101" spec="Taxon"/>
## <taxon id="sp102" spec="Taxon"/>

paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

## "sp15=2,sp24=2,sp25=1,sp27=1,sp29=2,sp34=2,sp36=2,sp37=2,sp39=2,sp40=2,sp41=2,sp42=1,sp43=2,sp46=2,sp47=2,sp48=1,sp49=2,sp50=2,sp51=2,sp54=1,sp55=1,sp56=2,sp57=2,sp59=2,sp60=2,sp61=2,sp62=2,sp63=2,sp68=2,sp69=2,sp70=2,sp71=1,sp72=2,sp73=2,sp76=2,sp77=2,sp78=2,sp79=1,sp81=2,sp82=2,sp83=2,sp84=2,sp85=2,sp86=2,sp87=2,sp88=1,sp89=2,sp90=2,sp91=2,sp92=2,sp93=2,sp94=2,sp95=2,sp96=2,sp97=2,sp98=2,sp99=2,sp100=2,sp101=2,sp102=2"

# finding mle
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars)

##      lambda0      lambda1          mu0          mu1          q01          q10 
## 2.304392e-08 3.622767e-01 4.535375e-01 6.484517e-02 2.127918e-01 1.526163e-01

## "(((sp93:0.2356622822,(sp94:0.2154739279,sp95:0.2154739279):0.02018835424):12.9847595,((((sp57:2.815062535,sp43:2.815062535):1.780249256,(sp36:4.449415911,((((sp99:0.06034332881,sp100:0.06034332881):0.4473268161,sp88:0.5076701449):1.726055497,sp56:2.233725641):1.371042764,((sp72:1.31183414,sp73:1.31183414):1.545878271,sp42:2.857712411):0.7470559947):0.8446475052):0.1458958804):3.207493326,(((sp27:5.033762494,(((((sp71:1.391331894,((sp86:0.5759627483,sp87:0.5759627483):0.7866153889,(sp89:0.4023765821,sp90:0.4023765821):0.9602015551):0.02875375666):1.595833121,(sp77:0.9006605643,sp78:0.9006605643):2.086504451):0.1392348633,sp41:3.126399879):0.8899363982,(((sp48:2.589801526,sp49:2.589801526):0.1590377796,sp46:2.748839306):0.6512535555,(sp60:2.03018895,((sp79:0.8625209504,(sp101:0.007232857647,sp102:0.007232857647):0.8552880928):0.6584553534,sp69:1.520976304):0.5092126459):1.369903912):0.6162434155):0.07444131832,sp34:4.090777595):0.9429848988):0.6018834209,sp25:5.635645915):0.003042611101,((sp59:2.160114586,(sp70:1.408818523,sp81:1.408818523):0.7512960633):2.006513375,(sp50:2.459826639,sp51:2.459826639):1.706801322):1.472060565):2.164116591):1.673667032,(((sp91:0.321558177,sp92:0.321558177):0.6062021719,sp76:0.9277603489):2.437759666,(sp47:2.651955001,(sp54:2.307131853,sp55:2.307131853):0.3448231482):0.7135650134):6.110952134):3.74394963):2.502755808,(sp15:15.70472892,(((sp37:3.930561202,(sp61:1.91883471,sp62:1.91883471):2.011726492):1.624114024,sp24:5.554675226):7.421080767,(((((sp84:0.6179056665,sp85:0.6179056665):0.2681327966,(sp82:0.8310677522,sp83:0.8310677522):0.05497071083):0.6598028498,sp68:1.545841313):3.050734938,sp29:4.596576251):2.433125762,(((sp96:0.1813892509,(sp97:0.1346470644,sp98:0.1346470644):0.04674218648):1.607626581,sp63:1.789015832):1.762305454,(sp39:3.456648649,sp40:3.456648649):0.0946726372):3.478380726):5.946053981):2.728972922):0.01844867133):0.0;"

# ----- BiSSE_fixed_tree_on_HiSSE_dataset.xml ----- #
pars <- c(.1,  .1,  .3,  # lambda 1, 2, 3
.05, .05, .05, # mu 1, 2, 3
.1, 0.0, # q12, q13
.1, .1, # q21, q23
0.0, .1 # q31, q32
) # pars above are equivalent to Fig. 1 in HiSSE paper

set.seed(10000)
phy <- tree.musse(pars, max.taxa=60, include.extinct=FALSE, x0=1)
phy$tip.state[phy$tip.state==3] <- 2 # hiding states
phy$tip.state <- phy$tip.state - 1

# tree
write.tree(phy)

## "((sp7:17.69971665,(sp34:2.921405715,sp35:2.921405715):14.77831094):15.46869175,(((((sp77:0.1860367355,sp78:0.1860367355):4.584657877,((sp67:0.8514863174,sp68:0.8514863174):0.8954959965,((sp64:1.077337234,(sp73:0.5016314604,sp74:0.5016314604):0.5757057734):0.1725445915,sp58:1.249881825):0.4971004886):3.023712299):7.768104514,sp25:12.53879913):6.812438184,((sp23:4.902320102,(sp55:1.527742454,((sp75:0.4873870334,sp76:0.4873870334):0.9538606252,(sp62:1.163335228,sp63:1.163335228):0.2779124303):0.08649479507):3.374577648):13.73148846,sp3:18.63380856):0.7174287518):7.688213403,(((((sp61:1.207186109,(sp71:0.6592772201,sp72:0.6592772201):0.5479088888):7.992858995,(((sp14:7.573004221,((sp31:3.946577942,((sp39:2.716827927,(sp59:1.223813007,sp60:1.223813007):1.49301492):0.2045427694,sp36:2.921370696):1.025207246):1.093651554,sp21:5.040229496):2.532774725):1.29161846,((((sp54:1.578469708,(sp69:0.6959006051,sp70:0.6959006051):0.8825691032):2.797784488,(sp56:1.512224659,sp57:1.512224659):2.864029537):0.5792204602,sp22:4.955474656):3.764087845,((sp51:2.056012337,(sp81:0.03182438438,sp82:0.03182438438):2.024187952):5.952203768,(sp38:3.212026203,sp33:3.212026203):4.796189902):0.7113463975):0.1450601792):0.0438553581,(((sp47:2.219906925,sp48:2.219906925):0.4858028937,sp40:2.705709819):4.257385058,((sp52:1.847137193,sp53:1.847137193):2.921544484,sp49:4.768681677):2.1944132):1.945383162):0.291567065):3.401238812,(sp41:2.700748896,(sp79:0.1442384287,sp80:0.1442384287):2.556510468):9.900535019):4.353908512,sp6:16.95519243):1.687331683,(((sp26:4.408240393,(sp45:2.255482777,sp46:2.255482777):2.152757616):2.8080253,sp15:7.216265693):5.259159804,((sp27:4.382339931,(sp43:2.385193162,sp44:2.385193162):1.997146769):3.208113605,((sp65:0.9805534352,sp66:0.9805534352):2.783753735,sp32:3.76430717):3.826146366):4.88497196):6.167098615):8.396926604):6.128957681):0.0;"

# taxa
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

## <taxon id="sp3" spec="Taxon"/>
## <taxon id="sp6" spec="Taxon"/>
## <taxon id="sp7" spec="Taxon"/>
## <taxon id="sp14" spec="Taxon"/>
## <taxon id="sp15" spec="Taxon"/>
## <taxon id="sp21" spec="Taxon"/>
## <taxon id="sp22" spec="Taxon"/>
## <taxon id="sp23" spec="Taxon"/>
## <taxon id="sp25" spec="Taxon"/>
## <taxon id="sp26" spec="Taxon"/>
## <taxon id="sp27" spec="Taxon"/>
## <taxon id="sp31" spec="Taxon"/>
## <taxon id="sp32" spec="Taxon"/>
## <taxon id="sp33" spec="Taxon"/>
## <taxon id="sp34" spec="Taxon"/>
## <taxon id="sp35" spec="Taxon"/>
## <taxon id="sp36" spec="Taxon"/>
## <taxon id="sp38" spec="Taxon"/>
## <taxon id="sp39" spec="Taxon"/>
## <taxon id="sp40" spec="Taxon"/>
## <taxon id="sp41" spec="Taxon"/>
## <taxon id="sp43" spec="Taxon"/>
## <taxon id="sp44" spec="Taxon"/>
## <taxon id="sp45" spec="Taxon"/>
## <taxon id="sp46" spec="Taxon"/>
## <taxon id="sp47" spec="Taxon"/>
## <taxon id="sp48" spec="Taxon"/>
## <taxon id="sp49" spec="Taxon"/>
## <taxon id="sp51" spec="Taxon"/>
## <taxon id="sp52" spec="Taxon"/>
## <taxon id="sp53" spec="Taxon"/>
## <taxon id="sp54" spec="Taxon"/>
## <taxon id="sp55" spec="Taxon"/>
## <taxon id="sp56" spec="Taxon"/>
## <taxon id="sp57" spec="Taxon"/>
## <taxon id="sp58" spec="Taxon"/>
## <taxon id="sp59" spec="Taxon"/>
## <taxon id="sp60" spec="Taxon"/>
## <taxon id="sp61" spec="Taxon"/>
## <taxon id="sp62" spec="Taxon"/>
## <taxon id="sp63" spec="Taxon"/>
## <taxon id="sp64" spec="Taxon"/>
## <taxon id="sp65" spec="Taxon"/>
## <taxon id="sp66" spec="Taxon"/>
## <taxon id="sp67" spec="Taxon"/>
## <taxon id="sp68" spec="Taxon"/>
## <taxon id="sp69" spec="Taxon"/>
## <taxon id="sp70" spec="Taxon"/>
## <taxon id="sp71" spec="Taxon"/>
## <taxon id="sp72" spec="Taxon"/>
## <taxon id="sp73" spec="Taxon"/>
## <taxon id="sp74" spec="Taxon"/>
## <taxon id="sp75" spec="Taxon"/>
## <taxon id="sp76" spec="Taxon"/>
## <taxon id="sp77" spec="Taxon"/>
## <taxon id="sp78" spec="Taxon"/>
## <taxon id="sp79" spec="Taxon"/>
## <taxon id="sp80" spec="Taxon"/>
## <taxon id="sp81" spec="Taxon"/>
## <taxon id="sp82" spec="Taxon"/>

# tip data
paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

## "sp3=2,sp6=2,sp7=1,sp14=2,sp15=1,sp21=2,sp22=2,sp23=2,sp25=2,sp26=1,sp27=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=2,sp36=2,sp38=2,sp39=2,sp40=2,sp41=2,sp43=2,sp44=2,sp45=2,sp46=1,sp47=2,sp48=2,sp49=2,sp51=2,sp52=2,sp53=2,sp54=2,sp55=2,sp56=2,sp57=2,sp58=2,sp59=2,sp60=2,sp61=2,sp62=2,sp63=2,sp64=2,sp65=2,sp66=2,sp67=2,sp68=2,sp69=2,sp70=2,sp71=2,sp72=2,sp73=2,sp74=2,sp75=2,sp76=2,sp77=2,sp78=2,sp79=2,sp80=2,sp81=1,sp82=1"

# finding MLE under BiSSE
pars.to.estimate <- pars <- c(.1, .3, .05, .05, .1, .1) # lambdas, mus, qs
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars)

##      lambda0      lambda1          mu0          mu1          q01          q10 
## 4.195862e-01 1.555696e-01 5.219573e-01 5.123616e-06 3.506821e-02 6.919963e-02 

# finding MLE under HiSSE
turnover.anc <- c(1,2,0,3)
eps.anc <- c(1,2,0,3)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates.nodual.no0B <- ParDrop(trans.rates,c(2,3,5,7,8,9,10,12))
trans.rates.nodual.no0B

sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)

pp <- hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.no0B, output.type="raw", root.type="equal", condition.on.survival=FALSE, root.p=NULL)

##      -lnL     AIC     AICc ntax
##  -182.084 384.168 388.6578   60

##                       lambda         mu
## rate0A             0.1594821 0.08712492
## rate1A             0.09142577 0.07480436
## rate0B                  0  0
## rate1B             0.2472443 0.1383931

## Transition Rates
##           (0A)       (1A) (0B)      (1B)
## (0A)        NA 0.29523936    0 0.0000000
## (1A) 0.4509993         NA    0 0.1355034
## (0B) 0.0000000 0.00000000   NA 0.0000000
## (1B) 0.0000000 0.01983649    0        NA
