# perform a simple analysis of sim study results

for (i in 1:3) {
	load(paste0("output/exp_",i,".RData"))

	N <- nrow(exp_res)
	s.good <-exp_res$s.success==TRUE
	ns.good<-exp_res$ns.success==TRUE
	good <- s.good&ns.good

	cat("==============\n")
	cat("Exp=",i,", Nreps=",N,", S good=",sum(s.good),", NS good=",sum(ns.good),", Both good=",sum(good),"\n",sep="")

	gres <- exp_res[good,]
	Ngood <- nrow(gres)

	cat("==============\n")
	cat("C LL:\n")
	cat("O = ",mean(gres$o.c_ll),", (", min(gres$o.c_ll), ", ", max(gres$o.c_ll), ")\n",sep="")
	cat("S = ",mean(gres$s.c_ll),", (", min(gres$s.c_ll), ", ", max(gres$s.c_ll), ")\n",sep="")
	cat("NS = ",mean(gres$ns.c_ll),", (", min(gres$ns.c_ll), ", ", max(gres$ns.c_ll), ")\n", sep="")
	cat("SE: O=",sd(gres$o.c_ll)/sqrt(Ngood),", S=",sd(gres$s.c_ll)/sqrt(Ngood),", NS=",sd(gres$ns.c_ll)/sqrt(Ngood),"\n",sep="")
	#print(t.test(gres$ns.c_ll,gres$s.c_ll,paired=T))

	cat("==============\n")
	cat("MSE:\n")
	cat("O = ",mean(gres$o.mse),", (", min(gres$o.mse), ", ", max(gres$o.mse), ")\n",sep="")
	cat("S = ",mean(gres$s.mse),", (", min(gres$s.mse), ", ", max(gres$s.mse), ")\n",sep="")
	cat("NS = ",mean(gres$ns.mse),", (", min(gres$ns.mse), ", ", max(gres$ns.mse), ")\n", sep="")
	cat("SE: O=",sd(gres$o.mse)/sqrt(Ngood),", S=",sd(gres$s.mse)/sqrt(Ngood),", NS=",sd(gres$ns.mse)/sqrt(Ngood),"\n",sep="")
	#print(t.test(gres$ns.mse,gres$s.mse,paired=T))
	#cat("NS MSE / S MSE = ",mean(gres$ns.mse)/mean(gres$s.mse),"\n",sep="")
	#print(t.test(gres$ns.mse/gres$s.mse))
	#cat("S MSE / NS MSE = ",mean(gres$s.mse)/mean(gres$ns.mse),"\n",sep="")
	#print(t.test(gres$s.mse/gres$ns.mse))

	cat("==============\n")
	cat("Computing time:\n")
	cat("O = ",mean(gres$o.elapsed),", (", min(gres$o.elapsed), ", ", max(gres$o.elapsed), ")\n",sep="")
	cat("S = ",mean(gres$s.elapsed),", (", min(gres$s.elapsed), ", ", max(gres$s.elapsed), ")\n",sep="")
	cat("NS = ",mean(gres$ns.elapsed)/60,", (", min(gres$ns.elapsed)/60, ", ", max(gres$ns.elapsed)/60, ")\n", sep="")
	cat("SE: O=",sd(gres$o.elapsed)/sqrt(Ngood),", S=",sd(gres$s.elapsed)/sqrt(Ngood),", NS=",sd(gres$ns.elapsed/60)/sqrt(Ngood),"\n",sep="")
	#print(t.test(gres$ns.elapsed/60,gres$s.elapsed/60,paired=T))

	cat("==============\n")
	cat("log(lambda) summary:\n")
	#print( summary(log(gres$ns.lambda)) )
	print( summary(gres$ns.lambda) )
}
