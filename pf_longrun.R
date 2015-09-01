# Example:  pair-formation script with extended (50 year) run
#
# Usage: R --vanilla < pf_longrun.R
#
# Rates are specified using monthly estimates
#
# Outputs are written to the output subdirectory.  They include:
#
#    lr_males.pdf - male population by disease state (long run)
#    lr_males_prev.pdf - same as a percentage of male population ("")
#    lr_females.pdf - female population by disease state ("")
#    lr_females_prev.pdf - same as a percentage of female population ("")
#    lr_pairs.pdf - paired population by disease state ("")
#    lr_pairs_prev.pdf - same as a percentage of pairs in population ("")
#    lr_overall.pdf - overall population by disease state ("")
#    lr_overall_prev.pdf - same as population percentages ("")
#
source("pf_model.R")
init<-init.pf(
              s.num=500000,
              s.num.g2=600000,
              i.num=37500,
              i.num.g2=60000,
              pairs.fract=0.4)
print(init)

param <-param.dcm(
              formation.rate=0.025,
              separate.rate=0.01,
              concur.rate=0.2,
              concur.rate.g2=0.05,
              b.rate=0.006,
              b.rate.g2=NA,
              ds.rate=0.0024,
              ds.rate.g2=0.0024,
              di.rate=0.012,
              di.rate.g2=0.008,
              trans.rate = 0.01,
              trans.rate.g2 = 0.04, 
	      balance="g1",
              p_act.rate = 4.0, 
              pi.rate = 0.25,
	      pfactor=5.,
              act.rate = 1.0)

control<-control.dcm(new.mod=pf_model,
	      type = "SI",
              nsteps = 600,
	      dt=0.1)

mod<-dcm(param,init,control)

#print(mod)

pdf("output/lr_males.pdf",width=8.5)
plot(mod, main="Male Population\n(50 year run)",y=c("s1.num","p1.num","i1.num","s2.num","p2.num","i2.num","num.g1"),popfrac=FALSE,
leg.name=c("single susc","single primary","single chronic","paired susc","paired primary","paired chronic","all males"),
col=c("blue","yellow","red","green","orange","pink","black"))
dev.off()

pdf("output/lr_males_prev.pdf",width=8.5)
plot(mod, main="Male Population (%)\n(50 year run)",y=c("prev_s.num.g1","prev_p.num.g1","prev_i.num.g1"),popfrac=FALSE,ylab="Percent",ylim=c(0,100),
leg.name=c("susceptible","primary stage","chronic stage"),
col=c("blue","yellow","red"))
dev.off()

pdf("output/lr_females.pdf",width=8.5)
plot(mod, main="Female Population\n(50 year run)",y=c("s1.num.g2","p1.num.g2","i1.num.g2","s2.num.g2","p2.num.g2","i2.num.g2","num.g2"),popfrac=FALSE,
leg.name=c("single susc","single primary","single chronic","paired susc","paired primary","paired chronic","all females"),
col=c("blue","yellow","red","green","orange","pink","black"))
dev.off()

pdf("output/lr_females_prev.pdf",width=8.5)
plot(mod, main="Female Population (%)\n(50 year run)",y=c("prev_s.num.g2","prev_p.num.g2","prev_i.num.g2"),popfrac=FALSE,ylab="Percent",ylim=c(0,100),
leg.name=c("susceptible","primary stage","chronic stage"),
col=c("blue","yellow","red"))
dev.off()

pdf("output/lr_pairs.pdf",width=8.5)
plot(mod, main="Paired Populations\n(50 year run)",y=c("p_ss.num","p_pp.num","p_sp.num","p_si.num","p_pi.num","p_ii.num","p_pairs.num"),popfrac=FALSE,
leg.name=c("both susc","both primary","susc / primary","susc / chronic","primary / chronic","both chronic","all pairs"),
col=c("blue","yellow","green","orange","pink","red","black"))
dev.off()

pdf("output/lr_pairs_prev.pdf",width=8.5)
plot(mod, main="Paired Populations (%)\n(50 year run)",y=c("prev_p_ss.num","prev_p_pp.num","prev_p_sp.num","prev_p_si.num","prev_p_pi.num","prev_p_ii.num"),popfrac=FALSE,ylab="Percent",ylim=c(0,100),
leg.name=c("both susc","both primary","susc / primary","susc / chronic","primary / chronic","both chronic"),
col=c("blue","yellow","green","orange","pink","red"))
dev.off()

pdf("output/lr_overall.pdf",width=8.5)
plot(mod, main="Total Population\n(50 year run)",y=c("t_s.num","t_p.num","t_i.num","t_num"),popfrac=FALSE,
leg.name=c("susceptible","primary","chronic","total"),
col=c("blue","green","red","black"))
dev.off()

pdf("output/lr_overall_prev.pdf",width=8.5)
plot(mod, main="Total Population (%)\n(50 year run)",y=c("prev_s.num","prev_p.num","prev_i.num"),popfrac=FALSE,
ylab="Percent",ylim=c(0,100),leg.name=c("susceptible","primary","chronic"),
col=c("blue","green","red"))
dev.off()
