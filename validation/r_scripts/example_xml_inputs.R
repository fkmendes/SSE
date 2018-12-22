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
pars <- c(.1,  .1,  .5,  # lambda 0, 1, 2
.05, .05, .05, # mu 0, 1, 2
.1, 0.0, # q01, q02
.1, .1, # q10, q12
0.0, .1 # q20, q21
) # pars above are equivalent to Fig. 1 in HiSSE paper

set.seed(10000)
phy <- tree.musse(pars, max.taxa=120, include.extinct=FALSE, x0=1)
phy$tip.state[phy$tip.state==3] <- 2 # hiding states
phy$tip.state <- phy$tip.state - 1

# tree
write.tree(phy)

## "((sp3:14.29757388,sp4:14.29757388)nd2:7.39121128,((((sp6:8.371777413,((((sp71:0.9344511607,sp72:0.9344511607)nd31:5.132807387,((sp24:3.034435126,((sp36:2.171713805,(sp95:0.4243839563,sp96:0.4243839563)nd93:1.747329849)nd92:0.008447696879,sp35:2.180161502)nd71:0.854273624)nd48:1.021856678,(sp61:1.285020052,sp62:1.285020052)nd49:2.771271752)nd32:2.010966743)nd18:1.179572333,(sp59:1.366915557,sp60:1.366915557)nd66:5.879915324)nd12:1.118832638,((((sp14:4.099159941,((sp132:0.05464629317,sp133:0.05464629317)nd95:3.919495433,(sp44:1.75101452,sp45:1.75101452)nd55:2.223127206)nd47:0.1250182145)nd40:0.5671972296,sp23:4.66635717)nd22:2.371860484,(sp41:2.03170906,((sp83:0.6340189991,(sp136:0.01502252686,sp137:0.01502252686)nd128:0.6189964722)nd122:0.3073797093,(sp130:0.09657770124,sp131:0.09657770124)nd123:0.8448210072)nd98:1.090310352)nd76:5.006508594)nd16:0.2990139288,((sp17:3.487998549,sp18:3.487998549)nd28:3.118610626,sp11:6.606609175)nd17:0.7306224088)nd13:1.028431935)nd11:0.006113894263)nd9:0.5530546196,(sp7:7.711799182,((((sp87:0.5852599838,sp88:0.5852599838)nd85:1.785935449,(sp122:0.1617929715,sp123:0.1617929715)nd86:2.209402461)nd37:3.070874501,(sp13:4.584734577,((((sp91:0.4850595759,sp92:0.4850595759)nd94:1.681337936,sp37:2.166397512)nd87:0.2039892456,sp30:2.370386757)nd69:0.7127652803,((sp134:0.02887418283,sp135:0.02887418283)nd117:1.038786103,(sp81:0.7003103016,sp82:0.7003103016)nd118:0.3673499838)nd70:2.015491752)nd42:1.501582539)nd38:0.8573353565)nd20:1.79370719,sp10:7.235777124)nd15:0.4760220587)nd10:1.21303285)nd6:3.460301764,(((((((sp25:2.968859494,((sp105:0.2791957709,(sp128:0.1044865585,sp129:0.1044865585)nd134:0.1747092123)nd91:1.928943011,sp34:2.208138782)nd72:0.7607207124)nd65:0.3833396916,sp20:3.352199186)nd63:0.03037964369,(sp106:0.2709886937,sp107:0.2709886937)nd64:3.111590136)nd45:0.904431931,(((((sp97:0.3846961595,sp98:0.3846961595)nd119:0.6765971808,(sp99:0.3766354031,sp100:0.3766354031)nd120:0.6846579372)nd79:1.382202644,(sp63:1.199303533,(sp67:1.105518883,(sp73:0.9159714603,((sp124:0.1354101187,sp125:0.1354101187)nd130:0.2080757167,sp103:0.3434858354)nd124:0.5724856249)nd114:0.1895474226)nd112:0.09378465058)nd80:1.244192451)nd75:0.5156567912,(sp64:1.157750038,sp65:1.157750038)nd74:1.801402737)nd67:0.2829626983,(sp120:0.1902877437,sp121:0.1902877437)nd68:3.05182773)nd46:1.044895286)nd33:1.779122799,(((sp32:2.210716154,sp33:2.210716154)nd89:1.816520108,((sp89:0.5024879985,sp90:0.5024879985)nd62:3.487660283,((((((sp93:0.4378159644,sp94:0.4378159644)nd121:0.5143800038,sp70:0.9521959682)nd108:0.4002410442,((sp110:0.244971814,sp111:0.244971814)nd113:0.8858630321,sp66:1.130834846)nd109:0.2216021663)nd105:0.1986776274,sp52:1.55111464)nd58:2.100451124,sp42:3.651565764)nd57:0.267867552,sp15:3.919433316)nd53:0.07071496549)nd51:0.03708798091)nd43:0.328603698,(sp28:2.453837079,sp29:2.453837079)nd44:1.902002882)nd34:1.710293599)nd29:0.4047538795,(((sp68:1.031949229,sp69:1.031949229)nd115:0.05637971444,(sp74:0.879501398,((sp112:0.2427729133,(sp116:0.2083130224,sp117:0.2083130224)nd135:0.03445989086)nd129:0.371298878,sp86:0.6140717913)nd125:0.2654296066)nd116:0.2088275454)nd81:1.334738238,((sp84:0.6279823274,sp85:0.6279823274)nd110:0.7239545343,(sp118:0.2067456007,sp119:0.2067456007)nd111:1.145191261)nd82:1.07113032)nd30:4.047820257)nd24:0.2231831415,((sp46:1.638872919,sp47:1.638872919)nd26:5.016081169,(((sp51:1.553137122,(sp55:1.420277503,((sp79:0.7135319078,sp80:0.7135319078)nd107:0.6626563229,sp58:1.376188231)nd106:0.04408927259)nd104:0.1328596191)nd60:2.015693403,(((((sp108:0.2585169656,sp109:0.2585169656)nd101:1.357595779,sp50:1.616112744)nd99:0.2651491555,(((sp101:0.3623348773,sp102:0.3623348773)nd127:0.3535188941,sp78:0.7158537714)nd102:0.8417642007,(sp77:0.7829898475,((((sp114:0.2341445922,sp115:0.2341445922)nd136:0.003211726445,sp113:0.2373563187)nd132:0.06707133001,(sp126:0.1324055725,sp127:0.1324055725)nd133:0.1720220762)nd131:0.006243845824,sp104:0.3106714945)nd126:0.472318353)nd103:0.7746281246)nd100:0.3236439277)nd83:0.5276629892,(sp75:0.8382681023,sp76:0.8382681023)nd84:1.570656787)nd77:0.06682142315,(sp39:2.073704257,(sp56:1.386573105,sp57:1.386573105)nd96:0.6871311522)nd78:0.4020420547)nd61:1.093084214)nd35:1.896776279,(sp31:2.347144917,sp48:2.347144917)nd36:3.118461888)nd27:1.189347283)nd25:0.03911649283)nd8:5.691063216)nd5:3.174693677,sp2:15.55982747)nd4:6.128957681):0.0;"

# taxa
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

## <taxon id="sp2" spec="Taxon"/>
## <taxon id="sp3" spec="Taxon"/>
## <taxon id="sp4" spec="Taxon"/>
## <taxon id="sp6" spec="Taxon"/>
## <taxon id="sp7" spec="Taxon"/>
## <taxon id="sp10" spec="Taxon"/>
## <taxon id="sp11" spec="Taxon"/>
## <taxon id="sp13" spec="Taxon"/>
## <taxon id="sp14" spec="Taxon"/>
## <taxon id="sp15" spec="Taxon"/>
## <taxon id="sp17" spec="Taxon"/>
## <taxon id="sp18" spec="Taxon"/>
## <taxon id="sp20" spec="Taxon"/>
## <taxon id="sp23" spec="Taxon"/>
## <taxon id="sp24" spec="Taxon"/>
## <taxon id="sp25" spec="Taxon"/>
## <taxon id="sp28" spec="Taxon"/>
## <taxon id="sp29" spec="Taxon"/>
## <taxon id="sp30" spec="Taxon"/>
## <taxon id="sp31" spec="Taxon"/>
## <taxon id="sp32" spec="Taxon"/>
## <taxon id="sp33" spec="Taxon"/>
## <taxon id="sp34" spec="Taxon"/>
## <taxon id="sp35" spec="Taxon"/>
## <taxon id="sp36" spec="Taxon"/>
## <taxon id="sp37" spec="Taxon"/>
## <taxon id="sp39" spec="Taxon"/>
## <taxon id="sp41" spec="Taxon"/>
## <taxon id="sp42" spec="Taxon"/>
## <taxon id="sp44" spec="Taxon"/>
## <taxon id="sp45" spec="Taxon"/>
## <taxon id="sp46" spec="Taxon"/>
## <taxon id="sp47" spec="Taxon"/>
## <taxon id="sp48" spec="Taxon"/>
## <taxon id="sp50" spec="Taxon"/>
## <taxon id="sp51" spec="Taxon"/>
## <taxon id="sp52" spec="Taxon"/>
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
## <taxon id="sp103" spec="Taxon"/>
## <taxon id="sp104" spec="Taxon"/>
## <taxon id="sp105" spec="Taxon"/>
## <taxon id="sp106" spec="Taxon"/>
## <taxon id="sp107" spec="Taxon"/>
## <taxon id="sp108" spec="Taxon"/>
## <taxon id="sp109" spec="Taxon"/>
## <taxon id="sp110" spec="Taxon"/>
## <taxon id="sp111" spec="Taxon"/>
## <taxon id="sp112" spec="Taxon"/>
## <taxon id="sp113" spec="Taxon"/>
## <taxon id="sp114" spec="Taxon"/>
## <taxon id="sp115" spec="Taxon"/>
## <taxon id="sp116" spec="Taxon"/>
## <taxon id="sp117" spec="Taxon"/>
## <taxon id="sp118" spec="Taxon"/>
## <taxon id="sp119" spec="Taxon"/>
## <taxon id="sp120" spec="Taxon"/>
## <taxon id="sp121" spec="Taxon"/>
## <taxon id="sp122" spec="Taxon"/>
## <taxon id="sp123" spec="Taxon"/>
## <taxon id="sp124" spec="Taxon"/>
## <taxon id="sp125" spec="Taxon"/>
## <taxon id="sp126" spec="Taxon"/>
## <taxon id="sp127" spec="Taxon"/>
## <taxon id="sp128" spec="Taxon"/>
## <taxon id="sp129" spec="Taxon"/>
## <taxon id="sp130" spec="Taxon"/>
## <taxon id="sp131" spec="Taxon"/>
## <taxon id="sp132" spec="Taxon"/>
## <taxon id="sp133" spec="Taxon"/>
## <taxon id="sp134" spec="Taxon"/>
## <taxon id="sp135" spec="Taxon"/>
## <taxon id="sp136" spec="Taxon"/>
## <taxon id="sp137" spec="Taxon"/>

# tip data
paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

## "sp2=1,sp3=2,sp4=1,sp6=2,sp7=2,sp10=2,sp11=2,sp13=1,sp14=2,sp15=2,sp17=1,sp18=1,sp20=2,sp23=2,sp24=2,sp25=2,sp28=2,sp29=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=2,sp36=2,sp37=2,sp39=2,sp41=2,sp42=2,sp44=2,sp45=2,sp46=2,sp47=2,sp48=2,sp50=2,sp51=2,sp52=2,sp55=2,sp56=2,sp57=2,sp58=2,sp59=2,sp60=2,sp61=2,sp62=2,sp63=2,sp64=2,sp65=2,sp66=2,sp67=2,sp68=2,sp69=2,sp70=2,sp71=2,sp72=2,sp73=2,sp74=2,sp75=2,sp76=2,sp77=2,sp78=2,sp79=2,sp80=2,sp81=2,sp82=2,sp83=2,sp84=2,sp85=2,sp86=2,sp87=2,sp88=2,sp89=2,sp90=2,sp91=2,sp92=2,sp93=2,sp94=2,sp95=2,sp96=2,sp97=2,sp98=2,sp99=2,sp100=2,sp101=2,sp102=2,sp103=2,sp104=2,sp105=2,sp106=2,sp107=2,sp108=2,sp109=2,sp110=2,sp111=2,sp112=2,sp113=2,sp114=2,sp115=2,sp116=2,sp117=2,sp118=2,sp119=2,sp120=2,sp121=2,sp122=2,sp123=2,sp124=2,sp125=2,sp126=2,sp127=2,sp128=2,sp129=2,sp130=2,sp131=2,sp132=2,sp133=2,sp134=2,sp135=2,sp136=2,sp137=2"

# finding MLE under BiSSE
pars.to.estimate <- pars <- c(.1, .5, .05, .05, .1, .1) # lambdas, mus, qs
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars)

##      lambda0      lambda1          mu0          mu1          q01          q10 
## 5.490635e-02 5.094458e-01 1.254943e-06 1.949697e-01 4.718035e-02 5.817885e-03 

# finding MLE under HiSSE
turnover.anc <- c(1,2,0,3)
eps.anc <- c(1,2,0,3)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates.nodual.no0B <- ParDrop(trans.rates,c(2,3,5,7,8,9,10,12))
trans.rates.nodual.no0B

sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)

pp <- hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.no0B, output.type="raw", root.type="equal", condition.on.survival=FALSE, root.p=NULL)

##      -lnL     AIC     AICc ntax
## -246.1055 512.211 514.2294  120

##                       lambda         mu
## rate0A             0.02913626 0.005376561
## rate1A             0.0636953 1.312858e-10
## rate1B             0.6093937 0.2787272

## Transition Rates
##            (0A)       (1A) (0B)       (1B)
## (0A)         NA 0.05197767    0 0.00000000
## (1A) 0.07369927         NA    0 0.03730514
## (0B) 0.00000000 0.00000000   NA 0.00000000
## (1B) 0.00000000 0.04861138    0         NA

# ----- BiSSE_fixed_tree_SCMfigure ----- #
## This is to make a figure examplifying SCM for the paper
pars <- c(.15, .3, .1, .1, .5, .5) # lambdas, mus, qs
set.seed(12345)
sampling.f<-c(1,1)
phy <- tree.bisse(pars, max.taxa=60, include.extinct=FALSE, x0=NA)
lik <- make.bisse(tree=phy, states=phy$tip.state, sampling.f=sampling.f, strict=FALSE)
anc.states <- asr.marginal(lik, pars)

## pal <- brewer.pal(8, "Set1")
## pal <- colorRampPalette(pal)(8)
## plot(phy, cex=.5, label.offset=0.2)
## nodelabels(pie=t(anc.states), cex=.5, piecol=c(pal[2],pal[6]))

write.tree(phy)

## "((sp29:18.53229859,(((sp53:3.131863819,sp60:3.131863819)nd38:10.94119267,(((((sp62:2.736899551,(sp73:1.874225183,(sp81:1.234108569,((sp91:0.808816581,(sp99:0.406087657,sp100:0.406087657)nd108:0.402728924)nd105:0.1366564909,sp86:0.9454730719)nd102:0.2886354976)nd92:0.6401166138)nd88:0.8626743675)nd83:0.1506048431,sp58:2.887504394)nd76:0.8497262898,(sp48:3.174400175,sp57:3.174400175)nd77:0.5628305084)nd69:2.831344081,(sp79:1.393112347,(sp82:1.178926953,(sp89:0.8338665622,sp90:0.8338665622)nd101:0.3450603904)nd96:0.2141853946)nd62:5.175462417)nd48:7.19199465,(sp41:8.826594412,(sp69:1.975074761,(sp70:1.956061282,(sp92:0.8075128508,sp93:0.8075128508)nd91:1.148548431)nd90:0.01901347888)nd57:6.851519651)nd35:4.933975002)nd33:0.3124870762)nd30:0.5213106024,(sp34:7.87385972,((sp102:0.3740821844,sp103:0.3740821844)nd66:5.227817761,(sp49:3.150298263,sp50:3.150298263)nd67:2.451601682)nd51:2.271959775)nd31:6.720507373)nd28:3.937931496)nd22:5.204457113,(((sp35:8.073085941,sp33:8.073085941)nd42:1.255836934,(((((sp94:0.787559482,(sp95:0.5938730482,sp96:0.5938730482)nd109:0.1936864338)nd107:0.06896461609,sp88:0.8565240981)nd93:1.952544229,sp59:2.809068327)nd78:0.4897100518,sp47:3.298778379)nd59:6.009321169,((sp68:2.163515355,(sp75:1.748864555,((sp97:0.4288077436,sp98:0.4288077436)nd103:0.5829152753,(sp106:0.2962868061,sp107:0.2962868061)nd104:0.7154362128)nd94:0.7371415358)nd89:0.4146508)nd65:5.622556397,(sp71:2.934221879,(sp104:0.2997879176,sp105:0.2997879176)nd82:2.634433962)nd53:4.851849872)nd45:1.522027797)nd43:0.02082332726)nd40:3.656198034,(((sp77:1.536978969,sp78:1.536978969)nd73:5.656794732,(((sp63:2.682789569,(sp76:1.740360696,((sp108:0.09727986203,sp109:0.09727986203)nd97:1.171891386,(((sp110:0.07842885097,sp111:0.07842885097)nd110:0.3134123023,sp101:0.3918411532)nd106:0.5508915421,sp87:0.9427326953)nd98:0.3264385523)nd95:0.4711894488)nd86:0.9424288723)nd74:1.091864496,(sp64:2.636433526,sp65:2.636433526)nd75:1.138220539)nd70:1.048060761,((sp84:1.002228203,sp85:1.002228203)nd99:0.2405533436,sp80:1.242781547)nd71:3.579933279)nd56:2.371058875)nd49:1.352103408,sp32:8.545877109)nd41:4.4392438)nd29:10.75163479)nd12:0.0;"

cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

## <taxon id="sp29" spec="Taxon"/>
## <taxon id="sp32" spec="Taxon"/>
## <taxon id="sp33" spec="Taxon"/>
## <taxon id="sp34" spec="Taxon"/>
## <taxon id="sp35" spec="Taxon"/>
## <taxon id="sp41" spec="Taxon"/>
## <taxon id="sp47" spec="Taxon"/>
## <taxon id="sp48" spec="Taxon"/>
## <taxon id="sp49" spec="Taxon"/>
## <taxon id="sp50" spec="Taxon"/>
## <taxon id="sp53" spec="Taxon"/>
## <taxon id="sp57" spec="Taxon"/>
## <taxon id="sp58" spec="Taxon"/>
## <taxon id="sp59" spec="Taxon"/>
## <taxon id="sp60" spec="Taxon"/>
## <taxon id="sp62" spec="Taxon"/>
## <taxon id="sp63" spec="Taxon"/>
## <taxon id="sp64" spec="Taxon"/>
## <taxon id="sp65" spec="Taxon"/>
## <taxon id="sp68" spec="Taxon"/>
## <taxon id="sp69" spec="Taxon"/>
## <taxon id="sp70" spec="Taxon"/>
## <taxon id="sp71" spec="Taxon"/>
## <taxon id="sp73" spec="Taxon"/>
## <taxon id="sp75" spec="Taxon"/>
## <taxon id="sp76" spec="Taxon"/>
## <taxon id="sp77" spec="Taxon"/>
## <taxon id="sp78" spec="Taxon"/>
## <taxon id="sp79" spec="Taxon"/>
## <taxon id="sp80" spec="Taxon"/>
## <taxon id="sp81" spec="Taxon"/>
## <taxon id="sp82" spec="Taxon"/>
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
## <taxon id="sp103" spec="Taxon"/>
## <taxon id="sp104" spec="Taxon"/>
## <taxon id="sp105" spec="Taxon"/>
## <taxon id="sp106" spec="Taxon"/>
## <taxon id="sp107" spec="Taxon"/>
## <taxon id="sp108" spec="Taxon"/>
## <taxon id="sp109" spec="Taxon"/>
## <taxon id="sp110" spec="Taxon"/>
## <taxon id="sp111" spec="Taxon"/>

paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

## "sp29=1,sp32=2,sp33=1,sp34=2,sp35=2,sp41=1,sp47=2,sp48=2,sp49=2,sp50=1,sp53=1,sp57=2,sp58=2,sp59=1,sp60=2,sp62=2,sp63=1,sp64=2,sp65=1,sp68=1,sp69=1,sp70=2,sp71=1,sp73=2,sp75=2,sp76=1,sp77=2,sp78=1,sp79=2,sp80=2,sp81=1,sp82=1,sp84=2,sp85=2,sp86=1,sp87=2,sp88=2,sp89=2,sp90=2,sp91=1,sp92=2,sp93=2,sp94=2,sp95=2,sp96=2,sp97=2,sp98=2,sp99=1,sp100=1,sp101=2,sp102=2,sp103=2,sp104=2,sp105=2,sp106=1,sp107=1,sp108=2,sp109=2,sp110=2,sp111=2"

# ----- FullMask_fixed_tree_on_HiSSE_RJHSDSEP_5x.xml ----- #
pars <- c(.1,  .1,  .5,  # lambda 1, 2, 3
.05, .05, .05, # mu 1, 2, 3
.1, 0.0, # q12, q13
.1, .1, # q21, q23
0.0, .1 # q31, q32
) # pars above are equivalent to Fig. 1 in HiSSE paper

set.seed(10000)
phy <- tree.musse(pars, max.taxa=120, include.extinct=FALSE, x0=1)
phy$tip.state[phy$tip.state==3] <- 2 # hiding states
phy$tip.state <- phy$tip.state - 1

# tree
write.tree(phy)

## "((sp3:14.29757388,sp4:14.29757388)nd2:7.39121128,((((sp6:8.371777413,((((sp71:0.9344511607,sp72:0.9344511607)nd31:5.132807387,((sp24:3.034435126,((sp36:2.171713805,(sp95:0.4243839563,sp96:0.4243839563)nd93:1.747329849)nd92:0.008447696879,sp35:2.180161502)nd71:0.854273624)nd48:1.021856678,(sp61:1.285020052,sp62:1.285020052)nd49:2.771271752)nd32:2.010966743)nd18:1.179572333,(sp59:1.366915557,sp60:1.366915557)nd66:5.879915324)nd12:1.118832638,((((sp14:4.099159941,((sp132:0.05464629317,sp133:0.05464629317)nd95:3.919495433,(sp44:1.75101452,sp45:1.75101452)nd55:2.223127206)nd47:0.1250182145)nd40:0.5671972296,sp23:4.66635717)nd22:2.371860484,(sp41:2.03170906,((sp83:0.6340189991,(sp136:0.01502252686,sp137:0.01502252686)nd128:0.6189964722)nd122:0.3073797093,(sp130:0.09657770124,sp131:0.09657770124)nd123:0.8448210072)nd98:1.090310352)nd76:5.006508594)nd16:0.2990139288,((sp17:3.487998549,sp18:3.487998549)nd28:3.118610626,sp11:6.606609175)nd17:0.7306224088)nd13:1.028431935)nd11:0.006113894263)nd9:0.5530546196,(sp7:7.711799182,((((sp87:0.5852599838,sp88:0.5852599838)nd85:1.785935449,(sp122:0.1617929715,sp123:0.1617929715)nd86:2.209402461)nd37:3.070874501,(sp13:4.584734577,((((sp91:0.4850595759,sp92:0.4850595759)nd94:1.681337936,sp37:2.166397512)nd87:0.2039892456,sp30:2.370386757)nd69:0.7127652803,((sp134:0.02887418283,sp135:0.02887418283)nd117:1.038786103,(sp81:0.7003103016,sp82:0.7003103016)nd118:0.3673499838)nd70:2.015491752)nd42:1.501582539)nd38:0.8573353565)nd20:1.79370719,sp10:7.235777124)nd15:0.4760220587)nd10:1.21303285)nd6:3.460301764,(((((((sp25:2.968859494,((sp105:0.2791957709,(sp128:0.1044865585,sp129:0.1044865585)nd134:0.1747092123)nd91:1.928943011,sp34:2.208138782)nd72:0.7607207124)nd65:0.3833396916,sp20:3.352199186)nd63:0.03037964369,(sp106:0.2709886937,sp107:0.2709886937)nd64:3.111590136)nd45:0.904431931,(((((sp97:0.3846961595,sp98:0.3846961595)nd119:0.6765971808,(sp99:0.3766354031,sp100:0.3766354031)nd120:0.6846579372)nd79:1.382202644,(sp63:1.199303533,(sp67:1.105518883,(sp73:0.9159714603,((sp124:0.1354101187,sp125:0.1354101187)nd130:0.2080757167,sp103:0.3434858354)nd124:0.5724856249)nd114:0.1895474226)nd112:0.09378465058)nd80:1.244192451)nd75:0.5156567912,(sp64:1.157750038,sp65:1.157750038)nd74:1.801402737)nd67:0.2829626983,(sp120:0.1902877437,sp121:0.1902877437)nd68:3.05182773)nd46:1.044895286)nd33:1.779122799,(((sp32:2.210716154,sp33:2.210716154)nd89:1.816520108,((sp89:0.5024879985,sp90:0.5024879985)nd62:3.487660283,((((((sp93:0.4378159644,sp94:0.4378159644)nd121:0.5143800038,sp70:0.9521959682)nd108:0.4002410442,((sp110:0.244971814,sp111:0.244971814)nd113:0.8858630321,sp66:1.130834846)nd109:0.2216021663)nd105:0.1986776274,sp52:1.55111464)nd58:2.100451124,sp42:3.651565764)nd57:0.267867552,sp15:3.919433316)nd53:0.07071496549)nd51:0.03708798091)nd43:0.328603698,(sp28:2.453837079,sp29:2.453837079)nd44:1.902002882)nd34:1.710293599)nd29:0.4047538795,(((sp68:1.031949229,sp69:1.031949229)nd115:0.05637971444,(sp74:0.879501398,((sp112:0.2427729133,(sp116:0.2083130224,sp117:0.2083130224)nd135:0.03445989086)nd129:0.371298878,sp86:0.6140717913)nd125:0.2654296066)nd116:0.2088275454)nd81:1.334738238,((sp84:0.6279823274,sp85:0.6279823274)nd110:0.7239545343,(sp118:0.2067456007,sp119:0.2067456007)nd111:1.145191261)nd82:1.07113032)nd30:4.047820257)nd24:0.2231831415,((sp46:1.638872919,sp47:1.638872919)nd26:5.016081169,(((sp51:1.553137122,(sp55:1.420277503,((sp79:0.7135319078,sp80:0.7135319078)nd107:0.6626563229,sp58:1.376188231)nd106:0.04408927259)nd104:0.1328596191)nd60:2.015693403,(((((sp108:0.2585169656,sp109:0.2585169656)nd101:1.357595779,sp50:1.616112744)nd99:0.2651491555,(((sp101:0.3623348773,sp102:0.3623348773)nd127:0.3535188941,sp78:0.7158537714)nd102:0.8417642007,(sp77:0.7829898475,((((sp114:0.2341445922,sp115:0.2341445922)nd136:0.003211726445,sp113:0.2373563187)nd132:0.06707133001,(sp126:0.1324055725,sp127:0.1324055725)nd133:0.1720220762)nd131:0.006243845824,sp104:0.3106714945)nd126:0.472318353)nd103:0.7746281246)nd100:0.3236439277)nd83:0.5276629892,(sp75:0.8382681023,sp76:0.8382681023)nd84:1.570656787)nd77:0.06682142315,(sp39:2.073704257,(sp56:1.386573105,sp57:1.386573105)nd96:0.6871311522)nd78:0.4020420547)nd61:1.093084214)nd35:1.896776279,(sp31:2.347144917,sp48:2.347144917)nd36:3.118461888)nd27:1.189347283)nd25:0.03911649283)nd8:5.691063216)nd5:3.174693677,sp2:15.55982747)nd4:6.128957681)nd1;"

cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

## <taxon id="sp2" spec="Taxon"/>
## <taxon id="sp3" spec="Taxon"/>
## <taxon id="sp4" spec="Taxon"/>
## <taxon id="sp6" spec="Taxon"/>
## <taxon id="sp7" spec="Taxon"/>
## <taxon id="sp10" spec="Taxon"/>
## <taxon id="sp11" spec="Taxon"/>
## <taxon id="sp13" spec="Taxon"/>
## <taxon id="sp14" spec="Taxon"/>
## <taxon id="sp15" spec="Taxon"/>
## <taxon id="sp17" spec="Taxon"/>
## <taxon id="sp18" spec="Taxon"/>
## <taxon id="sp20" spec="Taxon"/>
## <taxon id="sp23" spec="Taxon"/>
## <taxon id="sp24" spec="Taxon"/>
## <taxon id="sp25" spec="Taxon"/>
## <taxon id="sp28" spec="Taxon"/>
## <taxon id="sp29" spec="Taxon"/>
## <taxon id="sp30" spec="Taxon"/>
## <taxon id="sp31" spec="Taxon"/>
## <taxon id="sp32" spec="Taxon"/>
## <taxon id="sp33" spec="Taxon"/>
## <taxon id="sp34" spec="Taxon"/>
## <taxon id="sp35" spec="Taxon"/>
## <taxon id="sp36" spec="Taxon"/>
## <taxon id="sp37" spec="Taxon"/>
## <taxon id="sp39" spec="Taxon"/>
## <taxon id="sp41" spec="Taxon"/>
## <taxon id="sp42" spec="Taxon"/>
## <taxon id="sp44" spec="Taxon"/>
## <taxon id="sp45" spec="Taxon"/>
## <taxon id="sp46" spec="Taxon"/>
## <taxon id="sp47" spec="Taxon"/>
## <taxon id="sp48" spec="Taxon"/>
## <taxon id="sp50" spec="Taxon"/>
## <taxon id="sp51" spec="Taxon"/>
## <taxon id="sp52" spec="Taxon"/>
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
## <taxon id="sp103" spec="Taxon"/>
## <taxon id="sp104" spec="Taxon"/>
## <taxon id="sp105" spec="Taxon"/>
## <taxon id="sp106" spec="Taxon"/>
## <taxon id="sp107" spec="Taxon"/>
## <taxon id="sp108" spec="Taxon"/>
## <taxon id="sp109" spec="Taxon"/>
## <taxon id="sp110" spec="Taxon"/>
## <taxon id="sp111" spec="Taxon"/>
## <taxon id="sp112" spec="Taxon"/>
## <taxon id="sp113" spec="Taxon"/>
## <taxon id="sp114" spec="Taxon"/>
## <taxon id="sp115" spec="Taxon"/>
## <taxon id="sp116" spec="Taxon"/>
## <taxon id="sp117" spec="Taxon"/>
## <taxon id="sp118" spec="Taxon"/>
## <taxon id="sp119" spec="Taxon"/>
## <taxon id="sp120" spec="Taxon"/>
## <taxon id="sp121" spec="Taxon"/>
## <taxon id="sp122" spec="Taxon"/>
## <taxon id="sp123" spec="Taxon"/>
## <taxon id="sp124" spec="Taxon"/>
## <taxon id="sp125" spec="Taxon"/>
## <taxon id="sp126" spec="Taxon"/>
## <taxon id="sp127" spec="Taxon"/>
## <taxon id="sp128" spec="Taxon"/>
## <taxon id="sp129" spec="Taxon"/>
## <taxon id="sp130" spec="Taxon"/>
## <taxon id="sp131" spec="Taxon"/>
## <taxon id="sp132" spec="Taxon"/>
## <taxon id="sp133" spec="Taxon"/>
## <taxon id="sp134" spec="Taxon"/>
## <taxon id="sp135" spec="Taxon"/>
## <taxon id="sp136" spec="Taxon"/>
## <taxon id="sp137" spec="Taxon"/>

# tip data
paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")
