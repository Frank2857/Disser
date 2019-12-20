############################################################################
###################### MODELO 3 - PRIORI NÃO INFORMATIVA ###################
############################################################################
x = as.character(dados$ANO[3:23]);x
y= dados$FEV[3:23];y

cbind(x,y)

plot(x,y)

t0=12
n=length(y)
ncT=n-t0;ncT

############################################################################
# GEV COM TENDÊNCIA
sink("Gev_NOLINEAR_5_trend_NI.txt")
cat("
    model {

    # Verossimilhança
    for (i in 1:n) {
    y[i] ~ dgev(mi[i], sigma, eta)
    mi[i] <- mu-beta*exp(-gamma*x[i])
    }

    # Prior
    # dnorm(media, precisao)
    mu ~ dnorm(0,0.0001)
    sigma ~ dnorm(0,0.0001)
    eta ~ dunif(-1.5, 1.5)
    gamma ~ dunif(0.0, 1.0)
    beta  ~ dnorm(0.0,0.0001)

    #Para Fazer preditiva para Hipóses para a gamma H0=g<=g0
    gamma1<-0.1
    probH0g1<-step(gamma1-gamma)#Atribui valor cero, quando 'gamma1-gamma' menor ou aigual a cero '0.1'
    probH1g1<- 1-probH0g1       #Complementar de H0

    #Para Fazer preditiva para Hipóses para a beta H0=B<=B0
    beta0<-0
    probH0b<-step(beta0-beta) #Atribui valor cero, quando 'gamma1-gamma' menor ou aigual a cero '0'
    probH1b<- 1-probH0b       #Complementar de H0

    # Para Fazer preditiva
    #ncT<-9
    yp5 <- mu-beta*exp(-gamma*14) +((sigma/eta)*(pow(-log(1-1/5),-eta)-1))
    yp10 <- mu-beta*exp(-gamma*19) +((sigma/eta)*(pow(-log(1-1/10),-eta)-1))
    yp15 <- mu-beta*exp(-gamma*24) +((sigma/eta)*(pow(-log(1-1/15),-eta)-1))
    yp20 <- mu-beta*exp(-gamma*29) +((sigma/eta)*(pow(-log(1-1/20),-eta)-1))
    }
    ",fill=TRUE)
sink()

#plot(y)
trend=c(rep(0,t0),seq(1,ncT))
trend

dados_bug<- list(x=trend,y=y,n=length(y))
dados_bug

inits <- function(){ list(mu=33, sigma=10, beta=0.01, eta=0.01, gamma=0.01)}
params <- c("mu","sigma","eta","beta","gamma")
hipot<-c("probH0g1","probH1g1","probH0b","probH1b")
nc = 1      #Numero de cadeias
ni = 200000 #Tamanho da cadeira
nb = 50000  #Numero de simulação que serão descartadas
nt = 30     #salto (thin)

# Inicie o Amostrador
gev.bayes.trend_NL2 = bugs(data = dados_bug, inits = inits,
                           parameters =c(params,"yp5","yp10","yp15","yp20",hipot),
                           model = "Gev_NOLINEAR_5_trend_NI.txt",
                           n.thin = nt, n.chains = nc,
                           n.burnin = nb, n.iter = ni, codaPkg=FALSE, debug=T)

print(gev.bayes.trend_NL2, dig = 4)

post_gb_t_NL2<-as.mcmc(gev.bayes.trend_NL2$sims.matrix[,]) # salva a saída como cadeia mcmc

HPDinterval(post_gb_t_NL2) # Intervalo HPD

geweke.diag(post_gb_t_NL2)
raftery.diag(post_gb_t_NL2)
heidel.diag(post_gb_t_NL2)

par(mar=c(2,2,2,2))
plot(post_gb_t_NL2)


##############################################################
##############################################################
#       CALCULANDO MEDIDAS PARA DECISÃO PPARA BETA
##############################################################
##############################################################

pH0b=mean(post_gb_t_NL2[,12]);pH0b
pH1b=mean(post_gb_t_NL2[,13]);pH1b

########Evidencia
O_h1h0=pH1b/pH0b;O_h1h0# a chance de acontecer H1.


##############################################################
##############################################################
#       CALCULANDO MEDIDAS PARA DECISÃO PPARA GAMMA1=0.1
##############################################################
##############################################################

######Evidencia
pH0g1=mean(post_gb_t_NL2[,10]);pH0g1
pH1g1=mean(post_gb_t_NL2[,11]);pH1g1

O_h1h0=pH1g1/pH0g1;O_h1h0# a chance de acontecer H1

##############################################################
##############################################################
resumo2_NL2=print(gev.bayes.trend_NL2,dig=3) # salva a resumo da cadeia mcmc


#Preditiva
ypredl_t_NL2<-c(resumo2_NL2$mean$yp5,
                resumo2_NL2$mean$yp10,
                resumo2_NL2$mean$yp15,
                resumo2_NL2$mean$yp20) #salva as médias do y predito.
VPtNL2 = ypredl_t_NL2;VPtNL2
obs = c(34.6, 36.8, 36.8, 37.7)#Valores observados

EpGevT_NL2= abs((obs-VPtNL2)/obs)#Erro de predição
round(mean(EpGevT_NL2)*100,2)    #Erro médio de predição percentual.
