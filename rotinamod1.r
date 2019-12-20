############################################################################
###################### MODELO 1 - PRIORI NÃO INFORMATIVA ###################
############################################################################
x = as.character(dados$ANO[1:23]);x
y= dados$JAN[1:23];y
n=length(y)

############################################################################
# GEV COM TENDÊNCIA  LINEAR
sink("Gev_LINEAR_1_NI.txt")
cat("
    model {

    # Verossimilhança
    for (i in 1:n) {
    y[i] ~ dgev(mu, sigma, eta)
    }

    # Prior
    mu ~ dnorm(0,0.0001)
    sigma ~ dnorm(0,0.0001)
    eta ~ dunif(-1.5, 1.5)

    # Para Fazer preditiva
    yp5 <- mu + ((sigma/eta)*(pow(-log(1-1/5),-eta)-1))
    yp10 <- mu + ((sigma/eta)*(pow(-log(1-1/10),-eta)-1))
    yp15 <- mu + ((sigma/eta)*(pow(-log(1-1/15),-eta)-1))
    yp20 <- mu + ((sigma/eta)*(pow(-log(1-1/20),-eta)-1))
    }
    ",fill=TRUE)
sink()

dados_bug<- list(y=y,n=length(y))
dados_bug

inits <- function(){ list(mu=10, sigma=5, eta=0.01)}
params <- c("mu","sigma", "eta")
nc = 1      #Numero de cadeias
ni = 200000 #Tamanho da cadeira
nb = 50000  #Numero de simulação que serão descartadas
nt = 30     #Salto (thin)

# Inicie o Amostrador


gev.bayes_L1 = bugs(data = dados_bug, inits = inits,
                    parameters =c(params,"yp5","yp10","yp15","yp20"),
                    model = "Gev_LINEAR_1_NI.txt",
                    n.thin = nt, n.chains = nc,
                    n.burnin = nb, n.iter = ni, codaPkg=FALSE, debug=T)

print(gev.bayes_L1, dig = 4)

post_gb_L1<-as.mcmc(gev.bayes_L1$sims.matrix[,])#salva a saída como cadeia mcmc

geweke.diag(post_gb_L1)
raftery.diag(post_gb_L1)
heidel.diag(post_gb_L1)
par(mar=c(2,2,2,2))
plot(post_gb_L1)

HPDinterval(post_gb_L1)

resumo2_L1=print(gev.bayes_L1,dig=4) # salva a resumo da cadeia mcmc


#Preditiva
ypredl_L1<-c(resumo2_L1$mean$yp5,
             resumo2_L1$mean$yp10,
             resumo2_L1$mean$yp15,
             resumo2_L1$mean$yp20)#salva as médias do y predito.
VP_L1 = ypredl_L1;VP_L1
obs = c(36.4, 36.5, 36.5, 37.0)   #Valores observados

EpGev_L1= abs((obs-VP_L1)/obs)    #Erro de predição
round(mean(EpGev_L1)*100,2)       #Erro médio de predição percentual.
