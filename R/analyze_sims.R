# perform a simple analysis of sim study results

for (i in 3:3) {
	Nparts <- 20
	if (Nparts > 1) {
		cer <- data.frame()
		for (part in 1:Nparts) {
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

	cat("==============\n")
	cat("C LL:\n")
	cat("O = ",mean(gres$o.c_ll),", (", min(gres$o.c_ll), ", ", max(gres$o.c_ll), ")\n",sep="")
	cat("S = ",mean(gres$s.c_ll),", (", min(gres$s.c_ll), ", ", max(gres$s.c_ll), ")\n",sep="")
	cat("NS L1 = ",mean(gres$nsL1.c_ll),", (", min(gres$nsL1.c_ll), ", ", max(gres$nsL1.c_ll), ")\n", sep="")
	cat("NS L2 = ",mean(gres$nsL2.c_ll),", (", min(gres$nsL2.c_ll), ", ", max(gres$nsL2.c_ll), ")\n", sep="")
	cat("SE: O=",sd(gres$o.c_ll)/sqrt(Ngood),", S=",sd(gres$s.c_ll)/sqrt(Ngood),", NS L1=",sd(gres$nsL1.c_ll)/sqrt(Ngood),", NS L2=",sd(gres$nsL2.c_ll)/sqrt(Ngood),"\n",sep="")
	#print(t.test(gres$ns.c_ll,gres$s.c_ll,paired=T))

	cat("==============\n")
	cat("MSE:\n")
	cat("O = ",mean(gres$o.mse),", (", min(gres$o.mse), ", ", max(gres$o.mse), ")\n",sep="")
	cat("S = ",mean(gres$s.mse),", (", min(gres$s.mse), ", ", max(gres$s.mse), ")\n",sep="")
	cat("NS L1 = ",mean(gres$nsL1.mse),", (", min(gres$nsL1.mse), ", ", max(gres$nsL1.mse), ")\n", sep="")
	cat("NS L2 = ",mean(gres$nsL2.mse),", (", min(gres$nsL2.mse), ", ", max(gres$nsL2.mse), ")\n", sep="")
	cat("SE: O=",sd(gres$o.mse)/sqrt(Ngood),", S=",sd(gres$s.mse)/sqrt(Ngood),", NS L1=",sd(gres$nsL1.mse)/sqrt(Ngood),", NS L2=",sd(gres$nsL2.mse)/sqrt(Ngood),"\n",sep="")
	#print(t.test(gres$ns.mse,gres$s.mse,paired=T))
	#cat("NS MSE / S MSE = ",mean(gres$ns.mse)/mean(gres$s.mse),"\n",sep="")
	#print(t.test(gres$ns.mse/gres$s.mse))
	#cat("S MSE / NS MSE = ",mean(gres$s.mse)/mean(gres$ns.mse),"\n",sep="")
	#print(t.test(gres$s.mse/gres$ns.mse))

	cat("==============\n")
	cat("b0:\n")
	cat("O = ",mean(gres$o.b0),", (", min(gres$o.b0), ", ", max(gres$o.b0), ")\n",sep="")
	cat("S = ",mean(gres$s.b0),", (", min(gres$s.b0), ", ", max(gres$s.b0), ")\n",sep="")
	cat("NS L1 = ",mean(gres$nsL1.b0),", (", min(gres$nsL1.b0), ", ", max(gres$nsL1.b0), ")\n", sep="")
	cat("NS L2 = ",mean(gres$nsL2.b0),", (", min(gres$nsL2.b0), ", ", max(gres$nsL2.b0), ")\n", sep="")
	cat("SE: O=",sd(gres$o.b0)/sqrt(Ngood),", S=",sd(gres$s.b0)/sqrt(Ngood),", NS L1=",sd(gres$nsL1.b0)/sqrt(Ngood),", NS L2=",sd(gres$nsL2.b0)/sqrt(Ngood),"\n",sep="")

	cat("==============\n")
	cat("b1:\n")
	cat("O = ",mean(gres$o.b1),", (", min(gres$o.b1), ", ", max(gres$o.b1), ")\n",sep="")
	cat("S = ",mean(gres$s.b1),", (", min(gres$s.b1), ", ", max(gres$s.b1), ")\n",sep="")
	cat("NS L1 = ",mean(gres$nsL1.b1),", (", min(gres$nsL1.b1), ", ", max(gres$nsL1.b1), ")\n", sep="")
	cat("NS L2 = ",mean(gres$nsL2.b1),", (", min(gres$nsL2.b1), ", ", max(gres$nsL2.b1), ")\n", sep="")
	cat("SE: O=",sd(gres$o.b1)/sqrt(Ngood),", S=",sd(gres$s.b1)/sqrt(Ngood),", NS L1=",sd(gres$nsL1.b1)/sqrt(Ngood),", NS L2=",sd(gres$nsL2.b1)/sqrt(Ngood),"\n",sep="")

	cat("==============\n")
	cat("COV:\n")
	cat("O = ",mean(gres$o.cov),", (", min(gres$o.cov), ", ", max(gres$o.cov), ")\n",sep="")
	cat("S = ",mean(gres$s.cov),", (", min(gres$s.cov), ", ", max(gres$s.cov), ")\n",sep="")
	cat("NS L1 = ",mean(gres$nsL1.cov),", (", min(gres$nsL1.cov), ", ", max(gres$nsL1.cov), ")\n", sep="")
	cat("NS L2 = ",mean(gres$nsL2.cov),", (", min(gres$nsL2.cov), ", ", max(gres$nsL2.cov), ")\n", sep="")
	cat("SE: O=",sd(gres$o.cov)/sqrt(Ngood),", S=",sd(gres$s.cov)/sqrt(Ngood),", NS L1=",sd(gres$nsL1.cov)/sqrt(Ngood),", NS L2=",sd(gres$nsL2.cov)/sqrt(Ngood),"\n",sep="")

	cat("==============\n")
	cat("PI LEN:\n")
	cat("O = ",mean(gres$o.clen),", (", min(gres$o.clen), ", ", max(gres$o.clen), ")\n",sep="")
	cat("S = ",mean(gres$s.clen),", (", min(gres$s.clen), ", ", max(gres$s.clen), ")\n",sep="")
	cat("NS L1 = ",mean(gres$nsL1.clen),", (", min(gres$nsL1.clen), ", ", max(gres$nsL1.clen), ")\n", sep="")
	cat("NS L2 = ",mean(gres$nsL2.clen),", (", min(gres$nsL2.clen), ", ", max(gres$nsL2.clen), ")\n", sep="")
	cat("SE: O=",sd(gres$o.clen)/sqrt(Ngood),", S=",sd(gres$s.clen)/sqrt(Ngood),", NS L1=",sd(gres$nsL1.clen)/sqrt(Ngood),", NS L2=",sd(gres$nsL2.clen)/sqrt(Ngood),"\n",sep="")

	cat("==============\n")
	cat("Computing time:\n")
	cat("O = ",mean(gres$o.elapsed),", (", min(gres$o.elapsed), ", ", max(gres$o.elapsed), ")\n",sep="")
	cat("S = ",mean(gres$s.elapsed),", (", min(gres$s.elapsed), ", ", max(gres$s.elapsed), ")\n",sep="")
	cat("NS L1 = ",mean(gres$nsL1.elapsed)/60,", (", min(gres$nsL1.elapsed)/60, ", ", max(gres$nsL1.elapsed)/60, ")\n", sep="")
	cat("NS L2 = ",mean(gres$nsL2.elapsed)/60,", (", min(gres$nsL2.elapsed)/60, ", ", max(gres$nsL2.elapsed)/60, ")\n", sep="")
	cat("SE: O=",sd(gres$o.elapsed)/sqrt(Ngood),", S=",sd(gres$s.elapsed)/sqrt(Ngood),", NS L1=",sd(gres$nsL1.elapsed)/sqrt(Ngood),", NS L2=",sd(gres$nsL2.elapsed)/sqrt(Ngood),"\n",sep="")
	#print(t.test(gres$ns.elapsed/60,gres$s.elapsed/60,paired=T))

	cat("==============\n")
	cat("log(lambda) summary:\n")
	#print( summary(log(gres$ns.lambda)) )
	print( summary(gres$nsL1.lambda) )
	print( summary(gres$nsL2.lambda) )
}
