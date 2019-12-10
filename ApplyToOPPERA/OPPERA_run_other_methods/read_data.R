triggers.full <- read.csv("triggersnolydia.csv")

triggers <- triggers.full[,c(1,2,4,11,12,13,14,15,24,114,39:43,115,51,52,109,53,117,
                             110,111,60:64,119,81:86,125)]
	qhudate.r=rep(NA,nrow(triggers))
	enroldate.r=rep(NA,nrow(triggers))
	for (i in seq(nrow(triggers))){
	qhudate.r[i]=
	as.Date(toString(triggers$qhudatefilledout[i]), format="%d%b%Y")
	enroldate.r[i]=
	as.Date(toString(triggers$enrolmentdate[i]), format="%d%b%Y")
	timetoqhu=(qhudate.r-enroldate.r)/365.25
	}	
	triggers=cbind(triggers,qhudate.r,enroldate.r,timetoqhu)



triggers$siteid <- factor(triggers$siteid)
triggers$race <- factor(triggers$race)
triggers$gender <- factor(triggers$gender)
triggers$QHU_PainDuration <- factor(triggers$QHU_PainDuration)
triggers$QHU_NumMonthsPain <- factor(triggers$QHU_NumMonthsPain)
triggers$QHU_NumDaysPain <- factor(triggers$QHU_NumDaysPain)
