############################################################################
###################### MODELO 2 - PRIORI NÃO INFORMATIVA ###################
############################################################################
x = as.character(dados$ANO[1:23]);x
y= dados$FEV[1:23];y

t0=2          #Ponto do possível inicio de tendência
n=length(y)
ncT=n-t0;ncT

############################################################################
# GEV COM TENDÊNCIA  LINEAR
sink("Gev_TREND_LINEAR_NI.txt")
cat("
    model {

    # Verossimilhança
    for (i in 1:n) {
    y[i] ~ dgev(mi[i], sigma, eta)
    mi[i] <- mu + beta*x[i]
    }

    # Prior
    # dnorm(media, precisao)
    mu ~ dnorm(0,0.0001)
    sigma ~ dnorm(0,0.0001)
    eta ~ dunif(-1.5, 1.5)
    beta  ~ dunif(-5.0, 5.0)

    #Para Fazer preditiva para Hipóses para a beta H0=B<=B0
    beta0<-0
    probH0<-step(beta0-beta)# Atribui valor cero, quando 'beta0-beta' menor ou aigual a cero '0'
    probH1<- 1-probH0       # Complementar de H0
    # Para Fazer preditiva
    ncT<-21
    yp5 <- (mu + (beta*(ncT+5))) +((sigma/eta)*(pow(-log(1-1/5),-eta)-1))
    yp10 <- (mu + (beta*(ncT+10))) +((sigma/eta)*(pow(-log(1-1/10),-eta)-1))
    yp15 <- (mu + (beta*(ncT+15))) +((sigma/eta)*(pow(-log(1-1/15),-eta)-1))
    yp20 <- (mu + (beta*(ncT+20))) +((sigma/eta)*(pow(-log(1-1/20),-eta)-1))
    }
    ",fill=TRUE)
sink()

trend=c(rep(0,t0),seq(1,ncT))
trend

dados_bug<- list(x=trend,y=y,n=length(y))
dados_bug

inits <- function(){ list(mu=30, sigma=10, beta=0.01, eta=0.01)}
params <- c("mu","sigma","eta","beta")
nc = 1      #Numero de cadeias
ni = 200000 #Tamanho da cadeira
nb = 50000  #Numero de simulação que serão descartadas
nt = 30     #Salto (thin)

# Inicie o Amostrador
gev.bayes.trend_NL2 = bugs(data = dados_bug, inits = inits,
                           parameters =c(params,"probH0","probH1","yp5","yp10","yp15","yp20"),
                           model = "Gev_TREND_LINEAR_NI.txt",
                           n.thin = nt, n.chains = nc,
                           n.burnin = nb, n.iter = ni, codaPkg=FALSE, debug=T)

print(gev.bayes.trend_NL2, dig = 4)

post_gb_t_NL2<-as.mcmc(gev.bayes.trend_NL2$sims.matrix[,]) #salva a saída como cadeia mcmc

HPDinterval(post_gb_t_NL2) #Intervalo HPD

raftery.diag(post_gb_t_NL2)
geweke.diag(post_gb_t_NL2)
heidel.diag(post_gb_t_NL2)

par(mar=c(2,2,2,2))
plot(post_gb_t_NL2)

summary(post_gb_t_NL2 )

post_gb_t_NL2[1:20,] # imprimindo os 10 primeiros valores

#####################################################
#       CALCULANDO MEDIDAS PARA DECISÃO
##################################################
pH0=mean(post_gb_t_NL2[,5]);pH0
pH1=mean(post_gb_t_NL2[,6]);pH1

#  Evidencia
O_h1h0=pH1/pH0;O_h1h0 # A chance de acontecer H1

resumo2_NL2=print(gev.bayes.trend_NL2,dig=4) # salva a resumo da cadeia mcmc


#Preditiva
ypredl_t_NL2<-c(resumo2_NL2$mean$yp5,
                resumo2_NL2$mean$yp10,
                resumo2_NL2$mean$yp15,
                resumo2_NL2$mean$yp20) #salva as médias do y predito.
VPtNL2 = ypredl_t_NL2;VPtNL2
obs = c(34.6, 36.8, 36.8, 37.7)#Valores observados

EpGevT_NL2= abs((obs-VPtNL2)/obs) #Erro de predição
round(mean(EpGevT_NL2)*100,2)     #Erro médio de predição percentual.
