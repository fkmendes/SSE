library(ape)
library(diversitree)
#set.seed(1234)

args = commandArgs(TRUE)
exp.name = args[1]
num.states = as.integer(args[2])  # only accepts 4 or 8 right now

dir = "div/"
root = "/Users/jeff/Documents/Research/Phylogenetics/SSE/validation/"

# Set parameters
num.taxa = 22
num.trees = 100

symp.prob = 1.0
subsymp.prob = 1.0 / 6.0
vic.prob = 1.0 / 6.0
j.prob = 0

birth.rate = 0.4
death.rate = 0.1

s.rate = symp.prob * birth.rate
ss.rate = subsymp.prob * birth.rate
v.rate = vic.prob * birth.rate
j.rate = j.prob * birth.rate

classe.argnames = diversitree:::default.argnames.classe(num.states)
# "lambda111" "lambda112" "lambda113" "lambda114" "lambda122" "lambda123" "lambda124" "lambda133" "lambda134" "lambda144" "lambda211" "lambda212" "lambda213" "lambda214" "lambda222" "lambda223" "lambda224"
# "lambda233" "lambda234" "lambda244" "lambda311" "lambda312" "lambda313" "lambda314" "lambda322" "lambda323" "lambda324" "lambda333" "lambda334" "lambda344" "lambda411" "lambda412" "lambda413" "lambda414"
# "lambda422" "lambda423" "lambda424" "lambda433" "lambda434" "lambda444" "mu1"       "mu2"       "mu3"       "mu4"       "q12"       "q13"       "q14"       "q21"       "q23"       "q24"       "q31"    
# "q32"       "q34"       "q41"       "q42"       "q43"  

pars = numeric(length(classe.argnames))
names(pars) = classe.argnames
    
pars[names(pars) == "lambda111"] = s.rate
pars[names(pars) == "lambda222"] = s.rate
pars[names(pars) == "lambda333"] = s.rate
pars[names(pars) == "lambda223"] = j.rate
pars[names(pars) == "lambda323"] = j.rate
pars[names(pars) == "lambda423"] = v.rate
pars[names(pars) == "lambda424"] = ss.rate
pars[names(pars) == "lambda434"] = ss.rate

pars[names(pars) == "mu1"] = death.rate
pars[names(pars) == "mu2"] = death.rate
pars[names(pars) == "mu3"] = death.rate
pars[names(pars) == "mu4"] = death.rate

pars[names(pars) == "q12"] = 0.0
pars[names(pars) == "q13"] = 0.0
pars[names(pars) == "q14"] = 0.0
pars[names(pars) == "q21"] = 0.01
pars[names(pars) == "q23"] = 0.0
pars[names(pars) == "q24"] = 0.01
pars[names(pars) == "q31"] = 0.01
pars[names(pars) == "q32"] = 0.0
pars[names(pars) == "q34"] = 0.01
pars[names(pars) == "q41"] = 0.0
pars[names(pars) == "q42"] = 0.01
pars[names(pars) == "q43"] = 0.01


if (num.states == 8) {
    pars[names(pars) == "lambda555"] = s.rate
    pars[names(pars) == "lambda666"] = s.rate

    pars[names(pars) == "lambda556"] = j.rate
    pars[names(pars) == "lambda665"] = j.rate

    pars[names(pars) == "lambda765"] = v.rate
    pars[names(pars) == "lambda847"] = v.rate

    pars[names(pars) == "lambda775"] = ss.rate
    pars[names(pars) == "lambda776"] = ss.rate
    pars[names(pars) == "lambda884"] = ss.rate
    pars[names(pars) == "lambda887"] = ss.rate

    pars[names(pars) == "mu5"] = death.rate
    pars[names(pars) == "mu6"] = death.rate
    pars[names(pars) == "mu7"] = death.rate
    pars[names(pars) == "mu8"] = death.rate

    pars[names(pars) == "q51"] = 0.01
    pars[names(pars) == "q57"] = 0.01
    pars[names(pars) == "q61"] = 0.01
    pars[names(pars) == "q67"] = 0.01
    pars[names(pars) == "q75"] = 0.01
    pars[names(pars) == "q76"] = 0.01
    pars[names(pars) == "q84"] = 0.01
    pars[names(pars) == "q87"] = 0.01
}


# Get ground truth trees
start.state = num.states
phys = trees(pars, type="classe", max.taxa = num.taxa, x0=start.state, n=num.trees)

if (is.null(phys)) {
    print("bad tree")
    q()
}

ifelse(!dir.exists(dir), dir.create(dir), FALSE)




for (i in 1:length(phys)) {
    phy = phys[[i]]
    tips = phy$tip.state
    node.truth = phy$node.state
    
    
    print("Tree tips")
    print(tips)
    print("Tree in Newick format")
    write.tree(phy)
    print("Node truths")
    print(node.truth)
    
    
    tree.file.name = paste0(dir, exp.name, i, ".tree")
    write.tree(phy, file=tree.file.name)
    ntips = length(tips)
    
    
    node.truth.save = node.truth
    node.truth.names = names(node.truth)
    node.truth.save = t(node.truth.save)
    node.truth.save = rbind(node.truth.names, node.truth.save)
    node.file.name = paste0(dir, exp.name, i, "-node_truth.csv")
    write.table(node.truth.save, file=node.file.name, sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
    
    
    tips.save = tips
    tips.save.names = names(tips.save)
    tips.save = t(tips.save)
    tips.save = rbind(tips.save.names, tips.save)
    tips.file.name = paste0(dir, exp.name, i, "-tips.csv")
    write.table(tips.save, file=tips.file.name, sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
    
    
    beast.data = ""
    for (j in 1:length(tips)) {
        beast.data = paste0(beast.data, names(tips[j]), "=", tips[j], ",")
    }   
    beast.data = substr(beast.data, 0, nchar(beast.data) - 1)
    beast.data.file.name = paste0(root, dir, exp.name, i, "-beast_str.txt")
    cat(beast.data, file=beast.data.file.name)
    
    
    # Calculate likelyhood on tree
    if (num.states == 4) {
        sampling.f = c(1,1,1,1)
    } else if (num.states == 8) {
        sampling.f = c(1,1,1,1,1,1,1,1)
    }
    lik = make.classe(tree=phy, states=phy$tip.state, sampling.f=sampling.f, strict=FALSE, k=num.states)
    tree.lik = lik(pars=pars, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE) 
    # prior on root is the weighted average of D0 and D1, i.e., ROOT.OBS = D = D0 * (D0/(D0+D1)) + D1 * (D1/(D0+D1))
    print("Tree log likeyhood")
    print(tree.lik[1])
}


