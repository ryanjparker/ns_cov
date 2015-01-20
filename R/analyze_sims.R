# perform a simple analysis of sim study results

b0 <- 0
b1 <- 1

for (i in 3:3) {
	Nparts <- 20
	if (Nparts > 1) {
		cer <- data.frame()
		for (part in 1:Nparts) {
		#for (part in (1:Nparts)[-c(9,20)]) {
			load(paste0("output/exp_",i,"_",part,".RData"))
			cer <- rbind(cer, exp_res)
		}
		exp_res <- cer
	} else {
		load(paste0("output/exp_",i,".RData"))
	}

	N <- nrow(exp_res)
	s.good <-exp_res$s.success==TRUE
	nsL1.good<-exp_res$nsL1.success==TRUE
	nsL2.good<-exp_res$nsL2.success==TRUE
	good <- s.good&nsL1.good&nsL2.good

	cat("==============\n")
	cat("Exp=",i,", Nreps=",N,", S good=",sum(s.good),", NS L1 good=",sum(nsL1.good),", NS L2 good=",sum(nsL2.good),", Both good=",sum(good),"\n",sep="")

	gres <- exp_res[good,]
	Ngood <- nrow(gres)


# [1] "seed"           "n"              "o.elapsed"      "o.b0"          
# [5] "o.b1"           "o.b0.cov"       "o.b1.cov"       "o.b0.clen"     
# [9] "o.b1.clen"      "o.c_ll"         "o.mse"          "o.cov"         
#[13] "o.clen"         "s.success"      "s.elapsed"      "s.b0"          
#[17] "s.b1"           "s.b0.cov"       "s.b1.cov"       "s.b0.clen"     
#[21] "s.b1.clen"      "s.mse.tau"      "s.mse.sigma"    "s.mse.phi"     
#[25] "s.c_ll"         "s.mse"          "s.cov"          "s.clen"        
#[29] "nsL1.success"   "nsL1.elapsed"   "nsL1.b0"        "nsL1.b1"       
#[33] "nsL1.b0.cov"    "nsL1.b1.cov"    "nsL1.b0.clen"   "nsL1.b1.clen"  
#[37] "nsL1.mse.tau"   "nsL1.mse.sigma" "nsL1.mse.phi"   "nsL1.c_ll"     
#[41] "nsL1.mse"       "nsL1.cov"       "nsL1.clen"      "nsL1.lambda"   
#[45] "nsL2.success"   "nsL2.elapsed"   "nsL2.b0"        "nsL2.b1"       
#[49] "nsL2.b0.cov"    "nsL2.b1.cov"    "nsL2.b0.clen"   "nsL2.b1.clen"  
#[53] "nsL2.mse.tau"   "nsL2.mse.sigma" "nsL2.mse.phi"   "nsL2.c_ll"     
#[57] "nsL2.mse"       "nsL2.cov"       "nsL2.clen"      "nsL2.lambda"   

	gres_sum <- data.frame(t(colMeans(gres)))

	c_ll <- with(gres, list(s=t.test(o.c_ll,s.c_ll,paired=TRUE), nsL1=t.test(o.c_ll,nsL1.c_ll,paired=TRUE), nsL2=t.test(o.c_ll,nsL2.c_ll,paired=TRUE)))
	cat("==============\n")
	cat("C LL:\n")
	cat("O = ",with(gres_sum, round(o.c_ll-o.c_ll,2)),"\n",sep="")
	cat("S = ",with(gres_sum, round(o.c_ll-s.c_ll,2)),", (",with(c_ll$s, (conf.int[1] - estimate)/(qt(.05/2, df=parameter))),")\n",sep="")
	cat("NS L1 = ",with(gres_sum, round(o.c_ll-nsL1.c_ll,2)),", (",with(c_ll$nsL1, (conf.int[1] - estimate)/(qt(.05/2, df=parameter))),")\n",sep="")
	cat("NS L2 = ",with(gres_sum, round(o.c_ll-nsL2.c_ll,2)),", (",with(c_ll$nsL2, (conf.int[1] - estimate)/(qt(.05/2, df=parameter))),")\n",sep="")
	#cat("SE: O=",sd(gres$o.c_ll)/sqrt(Ngood),", S=",sd(gres$s.c_ll)/sqrt(Ngood),", NS L1=",sd(gres$nsL1.c_ll)/sqrt(Ngood),", NS L2=",sd(gres$nsL2.c_ll)/sqrt(Ngood),"\n",sep="")
	#print(t.test(gres$ns.c_ll,gres$s.c_ll,paired=T))

	cat("==============\n")
	cat("frob:\n")
	with(gres, cat("S = ",mean(s.frob),",",mean(s.frob)/mean(s.frob),"\n",sep=""))
	with(gres, cat("NS L1 = ",mean(nsL1.frob),", ",mean(nsL1.frob)/mean(s.frob),"\n",sep=""))
	with(gres, cat("NS L2 = ",mean(nsL2.frob),", ",mean(nsL2.frob)/mean(s.frob),"\n",sep=""))

	cat("==============\n")
	cat("rmse b0:\n")
	with(gres, cat("O = ",round(sqrt(mean((o.b0-b0)^2)),2),"\n",sep=""))
	with(gres, cat("S = ",round(sqrt(mean((s.b0-b0)^2)),2),", ",round(sqrt(mean((s.b0-b0)^2)/mean((o.b0-b0)^2)),2),"\n",sep=""))
	with(gres, cat("NS L1 = ",round(sqrt(mean((nsL1.b0-b0)^2)),2),", ",round(sqrt(mean((nsL1.b0-b0)^2)/mean((o.b0-b0)^2)),2),"\n",sep=""))
	with(gres, cat("NS L2 = ",round(sqrt(mean((nsL2.b0-b0)^2)),2),", ",round(sqrt(mean((nsL2.b0-b0)^2)/mean((o.b0-b0)^2)),2),"\n",sep=""))

	cat("==============\n")
	cat("rmse b1:\n")
	with(gres, cat("O = ",round(sqrt(mean((o.b1-b1)^2)),3),"\n",sep=""))
	with(gres, cat("S = ",round(sqrt(mean((s.b1-b1)^2)),3),", ",round(sqrt(mean((s.b1-b1)^2)/mean((o.b1-b1)^2)),2),"\n",sep=""))
	with(gres, cat("NS L1 = ",round(sqrt(mean((nsL1.b1-b1)^2)),3),", ",round(sqrt(mean((nsL1.b1-b1)^2)/mean((o.b1-b1)^2)),2),"\n",sep=""))
	with(gres, cat("NS L2 = ",round(sqrt(mean((nsL2.b1-b1)^2)),3),", ",round(sqrt(mean((nsL2.b1-b1)^2)/mean((o.b1-b1)^2)),2),"\n",sep=""))

	cat("==============\n")
	cat("cov b0:\n")
	with(gres, cat("O = ",mean(o.b0.cov),",\n",sep=""))
	with(gres, cat("S = ",mean(s.b0.cov),",\n",sep=""))
	with(gres, cat("NS L1 = ",mean(nsL1.b0.cov),",\n",sep=""))
	with(gres, cat("NS L2 = ",mean(nsL2.b0.cov),",\n",sep=""))

	cat("==============\n")
	cat("cov b1:\n")
	with(gres, cat("O = ",mean(o.b1.cov),",\n",sep=""))
	with(gres, cat("S = ",mean(s.b1.cov),",\n",sep=""))
	with(gres, cat("NS L1 = ",mean(nsL1.b1.cov),",\n",sep=""))
	with(gres, cat("NS L2 = ",mean(nsL2.b1.cov),",\n",sep=""))

	cat("==============\n")
	cat("irmse tau:\n")
	#with(gres, cat("O = ",round(sqrt(mean(o.mse.tau)),2),"\n",sep=""))
	with(gres, cat("S = ",round(sqrt(mean(s.mse.tau)),2),", ",round(sqrt(mean(s.mse.tau)/mean(s.mse.tau)),2),"\n",sep=""))
	with(gres, cat("NS L1 = ",round(sqrt(mean(nsL1.mse.tau)),2),", ",round(sqrt(mean(nsL1.mse.tau)/mean(s.mse.tau)),2),"\n",sep=""))
	with(gres, cat("NS L2 = ",round(sqrt(mean(nsL2.mse.tau)),2),", ",round(sqrt(mean(nsL2.mse.tau)/mean(s.mse.tau)),2),"\n",sep=""))

	cat("==============\n")
	cat("irmse sigma:\n")
	#with(gres, cat("O = ",round(sqrt(mean(o.mse.sigma)),2),"\n",sep=""))
	with(gres, cat("S = ",round(sqrt(mean(s.mse.sigma)),2),", ",round(sqrt(mean(s.mse.sigma)/mean(s.mse.sigma)),2),"\n",sep=""))
	with(gres, cat("NS L1 = ",round(sqrt(mean(nsL1.mse.sigma)),2),", ",round(sqrt(mean(nsL1.mse.sigma)/mean(s.mse.sigma)),2),"\n",sep=""))
	with(gres, cat("NS L2 = ",round(sqrt(mean(nsL2.mse.sigma)),2),", ",round(sqrt(mean(nsL2.mse.sigma)/mean(s.mse.sigma)),2),"\n",sep=""))

	cat("==============\n")
	cat("irmse phi:\n")
	#with(gres, cat("O = ",round(sqrt(mean(o.mse.phi)),2),"\n",sep=""))
	with(gres, cat("S = ",round(sqrt(mean(s.mse.phi)),3),", ",round(sqrt(mean(s.mse.phi)/mean(s.mse.phi)),2),"\n",sep=""))
	with(gres, cat("NS L1 = ",round(sqrt(mean(nsL1.mse.phi)),3),", ",round(sqrt(mean(nsL1.mse.phi)/mean(s.mse.phi)),2),"\n",sep=""))
	with(gres, cat("NS L2 = ",round(sqrt(mean(nsL2.mse.phi)),3),", ",round(sqrt(mean(nsL2.mse.phi)/mean(s.mse.phi)),2),"\n",sep=""))

	cat("==============\n")
	cat("rmse prediction:\n")
	with(gres, cat("O = ",round(sqrt(mean(o.mse)),3),"\n",sep=""))
	with(gres, cat("S = ",round(sqrt(mean(s.mse)),3),", ",round(sqrt(mean(s.mse)/mean(o.mse)),3),"\n",sep=""))
	with(gres, cat("NS L1 = ",round(sqrt(mean(nsL1.mse)),3),", ",round(sqrt(mean(nsL1.mse)/mean(o.mse)),3),"\n",sep=""))
	with(gres, cat("NS L2 = ",round(sqrt(mean(nsL2.mse)),3),", ",round(sqrt(mean(nsL2.mse)/mean(o.mse)),3),"\n",sep=""))

	cat("==============\n")
	cat("cov prediction:\n")
	with(gres, cat("O = ",mean(o.cov),",\n",sep=""))
	with(gres, cat("S = ",mean(s.cov),",\n",sep=""))
	with(gres, cat("NS L1 = ",mean(nsL1.cov),",\n",sep=""))
	with(gres, cat("NS L2 = ",mean(nsL2.cov),",\n",sep=""))

	cat("==============\n")
	cat("Computing time:\n")
	cat("O = ",mean(gres$o.elapsed),", (", min(gres$o.elapsed), ", ", max(gres$o.elapsed), ")\n",sep="")
	cat("S = ",mean(gres$s.elapsed),", (", min(gres$s.elapsed), ", ", max(gres$s.elapsed), ")\n",sep="")
	cat("NS L1 = ",mean(gres$nsL1.elapsed)/60,", (", min(gres$nsL1.elapsed)/60, ", ", max(gres$nsL1.elapsed)/60, ")\n", sep="")
	cat("NS L2 = ",mean(gres$nsL2.elapsed)/60,", (", min(gres$nsL2.elapsed)/60, ", ", max(gres$nsL2.elapsed)/60, ")\n", sep="")
	cat("SE: O=",sd(gres$o.elapsed)/sqrt(Ngood),", S=",sd(gres$s.elapsed)/sqrt(Ngood),", NS L1=",sd(gres$nsL1.elapsed)/sqrt(Ngood),", NS L2=",sd(gres$nsL2.elapsed)/sqrt(Ngood),"\n",sep="")
	#print(t.test(gres$ns.elapsed/60,gres$s.elapsed/60,paired=T))

	cat("==============\n")
	cat("lambda summary:\n")
	print( summary(gres$nsL1.lambda.tau) ); print( summary(gres$nsL2.lambda.tau) )
	print( summary(gres$nsL1.lambda.sigma) ); print( summary(gres$nsL2.lambda.sigma) )
	print( summary(gres$nsL1.lambda.phi) ); print( summary(gres$nsL2.lambda.phi) )
}
