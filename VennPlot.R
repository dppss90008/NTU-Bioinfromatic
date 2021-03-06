#### Root ####

H2O2 <- read.csv("H2O2.csv")


ABA <- read.csv("ABA_Shoot.csv")
GA <- read.csv("GA_Shoot.csv")
Auxin <- read.csv("Auxin_Shoot.csv")
BL <- read.csv("BL_Shoot.csv")
CK <- read.csv("CK_Shoot.csv")
JA <- read.csv("JA_Shoot.csv")


XD
venn(list("H2O2"=H2O2$name,"ABA"=ABA$name,"GA"=GA$name,"Auxin"=Auxin$name,"BL"=BL$name,"CK"=CK$name,"JA"=JA$name), ilabels = TRUE, 
    zcolor =c("dark blue","white","dark red","white","white","white","white"), size = 25, cexil = 1.2, cexsn = 1.5);

venn(list("H2O2"=H2O2$name,"ABA"=ABA$name,"JA"=JA$name,"CK"=CK$name,"Auxin"=Auxin$name), ilabels = TRUE, 
     zcolor = c("#003399","green","#F43E71","white","white"), size = 25, cexil = 1.2, cexsn = 1.5, opacity =  )

venn(list("H2O2"=H2O2$name,"GA"=GA$name,"BL"=BL$name), ilabels = TRUE,
     zcolor = "style", size = 25, cexil = 1.2, cexsn = 1.5,  )
area <- getZones("0-----1") # list of length 2
polygon(area[[1]], col="BLUE", opacity=0.3)

Intr <- intersect(H2O2$name,JA$name)
Intr <- intersect(Intr,ABA$name)
Intr <- intersect(Intr,Auxin$name)
Intr <- intersect(Intr,BL$name)
Intr <- intersect(Intr,CK$name)
Intr <- intersect(Intr,GA$name)

write.csv(Intr,file="H2O2_JA.csv")



# centroids for the two zones in the "E not A" zones
venn(5)
area <- getZones("0---1") # list of length 2
polygon(area[[1]], col="lightblue")
polygon(area[[2]], col="lightblue")
text(do.call("rbind", getCentroid(area)),
     labels = c("zone 1", "zone 2"), cex = 0.85)
