# ---------- #
# Trees and trait data for example .xmls
# ---------- #

library(ape)
library(diversitree)
library(hisse)

# ----- BiSSE_fixed_tree_SDSEP ----- #
pars <- c(.1, .5, .05, .05, .1, .1) # lambdas, mus, qs

set.seed(12345)
phy <- tree.bisse(pars, max.taxa=120, include.extinct=FALSE, x0=NA)

write.tree(phy)

## "((((((((sp88:0.7305723081,(sp110:0.3517943596,sp111:0.3517943596)nd133:0.3787779485)nd131:0.6126032115,((sp96:0.5196890137,sp97:0.5196890137)nd129:0.2622894529,sp78:0.7819784666)nd112:0.5611970529)nd75:0.8713987922,(sp57:1.455773138,((sp71:1.053741021,sp72:1.053741021)nd108:0.3540956393,((sp94:0.5321167068,sp95:0.5321167068)nd120:0.6052849744,(sp91:0.6118373316,sp92:0.6118373316)nd121:0.5255643496)nd109:0.2704349788)nd106:0.04793647802)nd76:0.7588011737)nd59:0.6454348926,(sp47:1.61994309,(sp50:1.532705744,(sp58:1.451997302,sp59:1.451997302)nd103:0.08070844236)nd98:0.08723734603)nd61:1.240066114)nd36:1.676965932,sp11:4.536975137)nd15:1.687950675,(((sp103:0.448272968,(sp138:0.08613692804,sp139:0.08613692804)nd138:0.36213604)nd137:0.009032491825,sp102:0.4573054599)nd35:4.399496701,sp6:4.856802161)nd26:1.368123651)nd5:1.152208972,((((sp33:2.24956865,sp35:2.24956865)nd46:1.426209219,((sp60:1.423855733,((sp104:0.4119903825,sp105:0.4119903825)nd139:0.02091855576,(sp134:0.1207488148,sp135:0.1207488148)nd140:0.3121601234)nd107:0.9909467946)nd96:0.2125091791,(((sp126:0.1645053988,sp127:0.1645053988)nd141:0.1835490384,sp112:0.3480544372)nd115:0.8698891337,(sp83:0.7454841432,sp84:0.7454841432)nd116:0.4724594277)nd97:0.418421341)nd66:2.039412957)nd9:3.023544015,(sp24:2.826445972,sp25:2.826445972)nd10:3.872875912)nd7:0.3077405585,((sp31:5.618597105,((sp119:0.2359467532,sp120:0.2359467532)nd82:1.904951537,sp46:2.14089829)nd28:3.477698814)nd11:0.7595809025,((sp53:1.485104655,sp54:1.485104655)nd93:0.2720998801,(sp115:0.285166619,sp116:0.285166619)nd94:1.472037916)nd12:4.620973472)nd8:0.6288844355)nd6:0.3700723412)nd2:0.8374129689,((sp42:1.88729009,sp43:1.88729009)nd16:4.224565776,((sp8:4.718334396,(((sp38:1.994098579,((sp85:0.7378233244,sp86:0.7378233244)nd118:0.4081087706,(sp87:0.7362545146,(sp130:0.1466536726,sp131:0.1466536726)nd132:0.5896008419)nd119:0.4096775804)nd87:0.8481664838)nd67:0.4437676967,((sp40:1.903126684,sp41:1.903126684)nd79:0.2787574296,sp34:2.181884113)nd68:0.2559821621)nd38:1.862291226,sp13:4.300157502)nd31:0.4181768941)nd30:0.8352288292,((((sp89:0.6569545471,(sp106:0.3787632156,sp107:0.3787632156)nd134:0.2781913315)nd48:2.975192865,sp19:3.632147412)nd45:0.1065481866,sp16:3.738695599)nd42:1.704978876,((((((sp140:0.05836369487,sp141:0.05836369487)nd110:1.315537159,sp63:1.373900854)nd91:0.4023286986,(sp44:1.712330625,(sp52:1.516571587,(sp113:0.3106816087,sp114:0.3106816087)nd105:1.205889978)nd95:0.1957590378)nd92:0.06389892751)nd77:0.4292478629,(sp66:1.177035351,sp67:1.177035351)nd78:1.028442064)nd51:1.221442003,(((sp22:2.858834248,sp23:2.858834248)nd57:0.162846485,((sp61:1.375699062,sp62:1.375699062)nd80:0.7656067357,(sp79:0.7762161959,sp80:0.7762161959)nd81:1.365089602)nd58:0.8803749348)nd55:0.1459512122,(sp69:1.06535337,sp70:1.06535337)nd56:2.102278575)nd52:0.2592874735)nd24:2.004495132,((sp142:0.02334854771,sp143:0.02334854771)nd39:4.244813798,(((sp27:2.611980347,(sp55:1.475452664,sp56:1.475452664)nd62:1.136527683)nd54:0.8918730549,((sp124:0.1853556377,sp125:0.1853556377)nd88:1.785377075,((sp48:1.561254771,(sp93:0.5806118125,(sp98:0.515688473,sp99:0.515688473)nd136:0.06492333946)nd101:0.9806429585)nd90:0.367401214,sp39:1.928655985)nd89:0.04207672764)nd50:1.53312069)nd43:0.2424696936,(sp20:3.244922908,((((sp51:1.51799645,((sp117:0.2761079111,sp118:0.2761079111)nd113:0.9590689502,((sp121:0.2176049823,sp122:0.2176049823)nd128:0.6643475616,sp77:0.8819525439)nd114:0.3532243174)nd104:0.2828195883)nd102:0.02918913232,sp49:1.547185582)nd73:0.7836372225,sp30:2.330822804)nd63:0.1182251958,(((sp65:2.106734091,(((sp75:0.9809616922,sp76:0.9809616922)nd122:0.01325772014,((sp132:0.1433480369,sp133:0.1433480369)nd142:0.0458249845,sp123:0.1891730214)nd123:0.8050463909)nd99:0.6120074427,((sp73:1.002219298,sp74:1.002219298)nd117:0.1558727881,sp68:1.158092086)nd100:0.4481347689)nd86:0.5005072356)nd71:0.2583449008,(sp37:2.112752022,((sp100:0.4954967354,sp101:0.4954967354)nd124:0.4618677482,((sp108:0.3550398,sp109:0.3550398)nd126:0.5564805588,((sp128:0.1587819134,sp129:0.1587819134)nd135:0.6047606251,sp81:0.7635425385)nd127:0.1479778203)nd125:0.04584412474)nd84:1.155387539)nd72:0.2523269692)nd69:0.01135267077,(sp136:0.09035357432,sp137:0.09035357432)nd70:2.286078088)nd64:0.07261633811)nd53:0.7958749077)nd44:0.5014001877)nd40:0.5218392502)nd25:1.163252205)nd23:0.01225992355)nd21:0.1098887506)nd17:0.5582926407)nd4:2.102691887):0.0;"

cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

## <taxon id="sp6" spec="Taxon"/>
## <taxon id="sp8" spec="Taxon"/>
## <taxon id="sp11" spec="Taxon"/>
## <taxon id="sp13" spec="Taxon"/>
## <taxon id="sp16" spec="Taxon"/>
## <taxon id="sp19" spec="Taxon"/>
## <taxon id="sp20" spec="Taxon"/>
## <taxon id="sp22" spec="Taxon"/>
## <taxon id="sp23" spec="Taxon"/>
## <taxon id="sp24" spec="Taxon"/>
## <taxon id="sp25" spec="Taxon"/>
## <taxon id="sp27" spec="Taxon"/>
## <taxon id="sp30" spec="Taxon"/>
## <taxon id="sp31" spec="Taxon"/>
## <taxon id="sp33" spec="Taxon"/>
## <taxon id="sp34" spec="Taxon"/>
## <taxon id="sp35" spec="Taxon"/>
## <taxon id="sp37" spec="Taxon"/>
## <taxon id="sp38" spec="Taxon"/>
## <taxon id="sp39" spec="Taxon"/>
## <taxon id="sp40" spec="Taxon"/>
## <taxon id="sp41" spec="Taxon"/>
## <taxon id="sp42" spec="Taxon"/>
## <taxon id="sp43" spec="Taxon"/>
## <taxon id="sp44" spec="Taxon"/>
## <taxon id="sp46" spec="Taxon"/>
## <taxon id="sp47" spec="Taxon"/>
## <taxon id="sp48" spec="Taxon"/>
## <taxon id="sp49" spec="Taxon"/>
## <taxon id="sp50" spec="Taxon"/>
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
## <taxon id="sp83" spec="Taxon"/>
## <taxon id="sp84" spec="Taxon"/>
## <taxon id="sp85" spec="Taxon"/>
## <taxon id="sp86" spec="Taxon"/>
## <taxon id="sp87" spec="Taxon"/>
## <taxon id="sp88" spec="Taxon"/>
## <taxon id="sp89" spec="Taxon"/>
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
## <taxon id="sp138" spec="Taxon"/>
## <taxon id="sp139" spec="Taxon"/>
## <taxon id="sp140" spec="Taxon"/>
## <taxon id="sp141" spec="Taxon"/>
## <taxon id="sp142" spec="Taxon"/>
## <taxon id="sp143" spec="Taxon"/>

paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

## "sp6=1,sp8=1,sp11=2,sp13=2,sp16=2,sp19=2,sp20=2,sp22=2,sp23=2,sp24=2,sp25=1,sp27=2,sp30=2,sp31=1,sp33=1,sp34=2,sp35=2,sp37=2,sp38=2,sp39=2,sp40=2,sp41=2,sp42=2,sp43=2,sp44=2,sp46=1,sp47=2,sp48=2,sp49=2,sp50=2,sp51=2,sp52=2,sp53=2,sp54=2,sp55=2,sp56=2,sp57=2,sp58=2,sp59=2,sp60=1,sp61=2,sp62=2,sp63=2,sp65=2,sp66=2,sp67=2,sp68=1,sp69=2,sp70=2,sp71=1,sp72=2,sp73=1,sp74=2,sp75=2,sp76=2,sp77=2,sp78=1,sp79=2,sp80=2,sp81=2,sp83=2,sp84=2,sp85=2,sp86=2,sp87=2,sp88=2,sp89=2,sp91=2,sp92=2,sp93=2,sp94=2,sp95=2,sp96=2,sp97=2,sp98=2,sp99=2,sp100=2,sp101=2,sp102=2,sp103=2,sp104=2,sp105=2,sp106=2,sp107=2,sp108=2,sp109=2,sp110=2,sp111=1,sp112=2,sp113=1,sp114=1,sp115=2,sp116=2,sp117=2,sp118=2,sp119=2,sp120=2,sp121=2,sp122=2,sp123=2,sp124=2,sp125=2,sp126=2,sp127=2,sp128=2,sp129=2,sp130=2,sp131=2,sp132=2,sp133=2,sp134=2,sp135=2,sp136=2,sp137=2,sp138=2,sp139=2,sp140=2,sp141=2,sp142=2,sp143=2"

# finding mle
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars)

##      lambda0      lambda1          mu0          mu1          q01          q10 
## 9.148731e-02 5.415631e-01 2.056091e-05 3.561561e-05 2.588266e-01 8.936368e-02 

# ----- ClaSSE_fixed_tree_SDSEP ----- #
argnames(make.classe(tree=phy, states=phy$tip.state+1, k=3, strict=FALSE)) # just to show order of parameters in classe object

##  [1] "lambda111" "lambda112" "lambda113" "lambda122" "lambda123" "lambda133"
##  [7] "lambda211" "lambda212" "lambda213" "lambda222" "lambda223" "lambda233"
## [13] "lambda311" "lambda312" "lambda313" "lambda322" "lambda323" "lambda333"
## [19] "mu1"       "mu2"       "mu3"       "q12"       "q13"       "q21"      
## [25] "q23"       "q31"       "q32"      

pars <- c(.1, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, .1, 0.0, 0.0,
          0.0, .5, .3, 0.0, .3, .1,
          .05, .05, .05, # mu1, mu2, mu3
          .1, .1, .1, # q12, q13, q21
          .1, .1, .1) # q23, q31, q32

set.seed(12345)
phy <- tree.classe(pars, max.taxa=120, x0=NA, include.extinct=FALSE)

# tree
write.tree(phy)

## "(((((((((sp151:0.2970300337,sp152:0.2970300337)nd161:0.416343853,sp142:0.7133738867)nd141:1.804695261,(sp92:2.407271845,sp93:2.407271845)nd142:0.110797303)nd90:3.136343939,sp41:5.654413087)nd53:2.295765957,((sp48:4.576332638,(sp50:4.480962459,((sp85:2.723247063,sp86:2.723247063)nd130:0.1364736545,(((sp145:0.5740892223,sp146:0.5740892223)nd153:0.5542578105,(sp129:1.102984711,sp130:1.102984711)nd154:0.02536232178)nd139:1.43897632,((sp101:2.207608979,sp102:2.207608979)nd145:0.188695966,sp94:2.396304945)nd140:0.1710184084)nd131:0.292397364)nd127:1.621241742)nd105:0.09537017956)nd55:3.344812061,sp28:7.921144699)nd54:0.0290343454)nd35:4.429668755,((sp117:1.460243154,sp118:1.460243154)nd108:5.913408404,(((sp126:1.176640433,sp127:1.176640433)nd109:3.27806327,((sp88:2.560187181,sp89:2.560187181)nd126:0.5030674068,sp71:3.063254588)nd110:1.391449115)nd63:2.675610383,(((sp155:0.22788305,sp156:0.22788305)nd163:0.1092294583,sp150:0.3371125083)nd158:0.6660442865,sp135:1.003156795)nd64:6.127157291)nd58:0.2433374723)nd36:5.006196242)nd34:0.1151394143,sp13:12.49498721)nd32:0.2806127433,(((sp36:6.07558979,(sp44:5.085275745,sp45:5.085275745)nd81:0.9903140443)nd47:3.555030305,sp20:9.630620094)nd44:0.8724231347,sp15:10.50304323)nd33:2.272556729)nd6:8.155879882,(((((((sp147:0.3650293315,(sp153:0.2845390462,sp154:0.2845390462)nd162:0.08049028529)nd147:1.305616776,(sp119:1.453329243,sp120:1.453329243)nd148:0.2173168652)nd124:1.452874128,sp69:3.123520236)nd99:3.328471099,(sp99:2.302653954,sp100:2.302653954)nd77:4.149337381)nd66:0.4252644714,((sp95:2.385849382,sp96:2.385849382)nd72:4.170162852,(((sp70:3.109966528,(sp82:2.741991085,sp83:2.741991085)nd125:0.3679754429)nd95:2.368968331,sp42:5.478934859)nd74:0.9882773095,((sp91:2.442187647,sp106:2.442187647)nd82:3.63321654,(sp157:0.2256171162,sp158:0.2256171162)nd83:5.849787071)nd80:0.3918079813)nd73:0.08880006541)nd67:0.3212435719)nd65:9.776438968,(((sp43:5.252015624,(sp79:2.803705655,(sp84:2.740623361,(sp113:1.562108161,(sp143:0.6598910829,sp144:0.6598910829)nd150:0.9022170779)nd135:1.1785152)nd132:0.06308229412)nd98:2.448309969)nd85:0.4993754086,(sp63:3.502885955,sp64:3.502885955)nd86:2.248505078)nd48:10.36466989,(sp51:4.46637739,sp52:4.46637739)nd106:11.64968353)nd18:0.5376338492)nd11:2.954773167,(((sp39:7.383926366,sp31:7.383926366)nd30:6.158442949,sp11:13.54236931)nd13:5.268821698,((sp18:16.58582574,(sp32:7.352791203,((sp121:1.434447698,(sp124:1.180741457,sp125:1.180741457)nd152:0.2537062412)nd146:5.377508553,(((sp138:0.8387875379,(sp140:0.8103231383,sp141:0.8103231383)nd159:0.02846439957)nd96:4.555656845,(sp90:2.497880285,((sp148:0.3590958059,sp149:0.3590958059)nd156:0.6825832171,(sp161:0.1209850427,sp162:0.1209850427)nd157:0.9206939803)nd143:1.456201262)nd97:2.896564098)nd70:1.269311365,(((sp115:1.488286205,sp116:1.488286205)nd89:4.18436034,sp40:5.672646545)nd87:0.04506981977,((sp133:1.031500626,sp134:1.031500626)nd133:1.750199466,((sp131:1.10123198,sp132:1.10123198)nd136:1.572715306,(sp128:1.127079621,(sp136:0.9794983809,sp137:0.9794983809)nd155:0.1475812398)nd149:1.546867665)nd134:0.1077528064)nd119:2.936016273)nd71:0.9460393827)nd69:0.1482005039)nd59:0.5408349518)nd22:9.233034535)nd15:1.849852712,(((((sp57:3.920283848,((sp107:1.813907675,sp108:1.813907675)nd117:1.913958592,(sp139:0.8344412231,(sp163:0.0357731472,sp164:0.0357731472)nd160:0.7986680759)nd123:2.893425044)nd114:0.1924175807)nd101:1.0115441,((sp60:3.578325358,(sp103:2.115068739,sp104:2.115068739)nd138:1.463256619)nd103:1.340745351,sp65:4.919070709)nd102:0.01275723948)nd49:6.19690397,(((sp97:5.497329113,(sp114:1.513451635,(sp122:1.362934701,sp123:1.362934701)nd151:0.1505169334)nd94:3.983877479)nd60:1.738287931,(((sp73:2.996910825,sp78:2.996910825)nd115:0.7421721701,(sp80:2.752739442,sp81:2.752739442)nd116:0.9863435536)nd91:1.813672423,(sp110:1.731646589,sp111:1.731646589)nd111:3.821108829)nd84:1.682861625)nd42:3.474648558,(sp46:4.714759442,sp47:4.714759442)nd51:5.99550616)nd41:0.4184663164)nd31:3.858130051,(((sp56:3.943140603,(sp159:0.1212999078,sp160:0.1212999078)nd113:3.821840695)nd78:2.192987178,(sp54:4.092328803,sp55:4.092328803)nd79:2.043798978)nd50:4.672985872,(sp29:7.566981444,sp30:7.566981444)nd40:3.242132208)nd29:4.177748316)nd25:0.9196375993,(sp26:8.607195636,sp27:8.607195636)nd27:7.299303933)nd16:2.529178882)nd14:0.3755125623)nd12:0.797276928)nd10:1.323011899):0.0;"

# taxa
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

## <taxon id="sp11" spec="Taxon"/>
## <taxon id="sp13" spec="Taxon"/>
## <taxon id="sp15" spec="Taxon"/>
## <taxon id="sp18" spec="Taxon"/>
## <taxon id="sp20" spec="Taxon"/>
## <taxon id="sp26" spec="Taxon"/>
## <taxon id="sp27" spec="Taxon"/>
## <taxon id="sp28" spec="Taxon"/>
## <taxon id="sp29" spec="Taxon"/>
## <taxon id="sp30" spec="Taxon"/>
## <taxon id="sp31" spec="Taxon"/>
## <taxon id="sp32" spec="Taxon"/>
## <taxon id="sp36" spec="Taxon"/>
## <taxon id="sp39" spec="Taxon"/>
## <taxon id="sp40" spec="Taxon"/>
## <taxon id="sp41" spec="Taxon"/>
## <taxon id="sp42" spec="Taxon"/>
## <taxon id="sp43" spec="Taxon"/>
## <taxon id="sp44" spec="Taxon"/>
## <taxon id="sp45" spec="Taxon"/>
## <taxon id="sp46" spec="Taxon"/>
## <taxon id="sp47" spec="Taxon"/>
## <taxon id="sp48" spec="Taxon"/>
## <taxon id="sp50" spec="Taxon"/>
## <taxon id="sp51" spec="Taxon"/>
## <taxon id="sp52" spec="Taxon"/>
## <taxon id="sp54" spec="Taxon"/>
## <taxon id="sp55" spec="Taxon"/>
## <taxon id="sp56" spec="Taxon"/>
## <taxon id="sp57" spec="Taxon"/>
## <taxon id="sp60" spec="Taxon"/>
## <taxon id="sp63" spec="Taxon"/>
## <taxon id="sp64" spec="Taxon"/>
## <taxon id="sp65" spec="Taxon"/>
## <taxon id="sp69" spec="Taxon"/>
## <taxon id="sp70" spec="Taxon"/>
## <taxon id="sp71" spec="Taxon"/>
## <taxon id="sp73" spec="Taxon"/>
## <taxon id="sp78" spec="Taxon"/>
## <taxon id="sp79" spec="Taxon"/>
## <taxon id="sp80" spec="Taxon"/>
## <taxon id="sp81" spec="Taxon"/>
## <taxon id="sp82" spec="Taxon"/>
## <taxon id="sp83" spec="Taxon"/>
## <taxon id="sp84" spec="Taxon"/>
## <taxon id="sp85" spec="Taxon"/>
## <taxon id="sp86" spec="Taxon"/>
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
## <taxon id="sp99" spec="Taxon"/>
## <taxon id="sp100" spec="Taxon"/>
## <taxon id="sp101" spec="Taxon"/>
## <taxon id="sp102" spec="Taxon"/>
## <taxon id="sp103" spec="Taxon"/>
## <taxon id="sp104" spec="Taxon"/>
## <taxon id="sp106" spec="Taxon"/>
## <taxon id="sp107" spec="Taxon"/>
## <taxon id="sp108" spec="Taxon"/>
## <taxon id="sp110" spec="Taxon"/>
## <taxon id="sp111" spec="Taxon"/>
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
## <taxon id="sp138" spec="Taxon"/>
## <taxon id="sp139" spec="Taxon"/>
## <taxon id="sp140" spec="Taxon"/>
## <taxon id="sp141" spec="Taxon"/>
## <taxon id="sp142" spec="Taxon"/>
## <taxon id="sp143" spec="Taxon"/>
## <taxon id="sp144" spec="Taxon"/>
## <taxon id="sp145" spec="Taxon"/>
## <taxon id="sp146" spec="Taxon"/>
## <taxon id="sp147" spec="Taxon"/>
## <taxon id="sp148" spec="Taxon"/>
## <taxon id="sp149" spec="Taxon"/>
## <taxon id="sp150" spec="Taxon"/>
## <taxon id="sp151" spec="Taxon"/>
## <taxon id="sp152" spec="Taxon"/>
## <taxon id="sp153" spec="Taxon"/>
## <taxon id="sp154" spec="Taxon"/>
## <taxon id="sp155" spec="Taxon"/>
## <taxon id="sp156" spec="Taxon"/>
## <taxon id="sp157" spec="Taxon"/>
## <taxon id="sp158" spec="Taxon"/>
## <taxon id="sp159" spec="Taxon"/>
## <taxon id="sp160" spec="Taxon"/>
## <taxon id="sp161" spec="Taxon"/>
## <taxon id="sp162" spec="Taxon"/>
## <taxon id="sp163" spec="Taxon"/>
## <taxon id="sp164" spec="Taxon"/>

# tip data
paste(paste(phy$tip.label, (phy$tip.state), sep="="), collapse=",")

## "sp11=2,sp13=1,sp15=1,sp18=1,sp20=3,sp26=3,sp27=2,sp28=2,sp29=1,sp30=2,sp31=2,sp32=2,sp36=2,sp39=3,sp40=2,sp41=1,sp42=1,sp43=2,sp44=2,sp45=2,sp46=3,sp47=1,sp48=3,sp50=2,sp51=1,sp52=2,sp54=2,sp55=2,sp56=2,sp57=1,sp60=2,sp63=1,sp64=3,sp65=2,sp69=3,sp70=1,sp71=2,sp73=2,sp78=3,sp79=2,sp80=1,sp81=1,sp82=1,sp83=2,sp84=1,sp85=1,sp86=2,sp88=3,sp89=1,sp90=1,sp91=2,sp92=1,sp93=2,sp94=3,sp95=1,sp96=1,sp97=3,sp99=1,sp100=2,sp101=1,sp102=1,sp103=2,sp104=2,sp106=2,sp107=1,sp108=1,sp110=2,sp111=2,sp113=2,sp114=2,sp115=2,sp116=2,sp117=2,sp118=2,sp119=1,sp120=2,sp121=1,sp122=1,sp123=2,sp124=1,sp125=1,sp126=2,sp127=3,sp128=2,sp129=1,sp130=2,sp131=1,sp132=1,sp133=1,sp134=1,sp135=2,sp136=3,sp137=3,sp138=2,sp139=1,sp140=1,sp141=2,sp142=3,sp143=2,sp144=3,sp145=1,sp146=3,sp147=2,sp148=1,sp149=1,sp150=1,sp151=1,sp152=2,sp153=1,sp154=2,sp155=1,sp156=1,sp157=1,sp158=2,sp159=2,sp160=2,sp161=2,sp162=2,sp163=1,sp164=3"

# finding mle
lik <- make.classe(phy, phy$tip.state, k=3)
fit <- find.mle(lik, pars) # bombs!

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
