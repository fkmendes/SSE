library(RColorBrewer,ggpubr)

p.a <- get_palette(c("#FFFFFF", "#FF00DC"), k=5)
p.a
# #FFFFFF #FFBFF6 #FF7FED #FF3FE4 #FF00DC
# 0,0,0,0 0,3,28,0, 11,56,0,0 18,78,0,0 18,85,0,0

p.c <- get_palette(c("#FFFFFF", "#F1C40F"), k=5)
p.c
# #FFFFFF #FBF0C3 #F8E187 #F4D24B #F1C40F
# 0,0,0,0 0,4,22,2 0,9,44,3 0,13,66,4 0,18,89,5

p.g <- get_palette(c("#FFFFFF", "#34495E"), k=5)
p.g
# #FFFFFF #CCD1D6 #99A3AE #667686 #34495E
# 0,0,0,0 4,2,0,16 8,4,0,32 13,6,0,47 16,8,0,63

p.t <- get_palette(c("#FFFFFF", "#3498DB"), k=5)
p.t
# #FFFFFF #CCE5F6 #99CBED #66B1E3 #3498DB
# 0,0,0,0 16,7,0,4 33,13,0,7 49,20,0,11 65,26,0,14

p.q <- get_palette(c("#FFFFFF", "#C2C1C0"), k=5)
p.q
# #FFFFFF #EFEFEF #E0E0DF #D1D0CF #C2C1C0
# 0,0,0,0 0,0,0,6 0,0,0,12 0,0,1,18 0,0,1,24
