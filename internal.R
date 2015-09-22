


GetDataFromBIOM <-function(dataBIOM)
{
  
  counts = biom_data(dataBIOM)
  taxo = observation_metadata(dataBIOM)
  return(list(counts=counts,taxo=taxo))
}


GetDataFromCT <-function(dataC,dataT)
{
  
  counts = dataC
  taxo = dataT
  return(list(counts=counts,taxo=taxo))
}




DescriptiveStat <-function(vect)
{
  
  nbmiss = length(which(is.na((vect))))
  nbval = length(vect) - nbmiss
  sum = sum(vect,na.rm=TRUE)
  moy = mean(vect,na.rm=TRUE)
  var = var(vect,na.rm=TRUE)
  sd = sd(vect,na.rm=TRUE)
  CV = sd/moy
  stat = c(nbval,nbmiss,sum, summary(vect),var,sd,CV)
  
  names(stat) = c("Nb valeurs","Nb manquants","Somme","Min",
                  "1er Quartile","Mediane",'Moyenne',"3eme Quartile","Max","Variance","Ecart-type","Coeff Variation")
    
  return(stat)
}

DescriptiveStatQuali <-function(dataQuali,namesQuali,indic)
{
  
  res = matrix(0,ncol=5)
  
  for(i in 1:ncol(dataQuali))
  {
    datatmp = dataQuali[,i]
    tabQuali = table(datatmp)
    res2 = c(names(tabQuali),"Total")
    res1 = c(namesQuali[i],rep(NA,length(res2)-1))
    res3 = c(tabQuali,sum(tabQuali))
    res4 = signif(c(tabQuali/sum(tabQuali),1)*100,2)
    res5 = c(signif(cumsum(tabQuali/sum(tabQuali))*100,2),NA)
    res = rbind(res,cbind(res1,res2,res3,res4,res5))
  }
  res=res[-1,]
  rownames(res)=res[,1]
  res=res[,-1]
  colnames(res) = c("Modalites","Effectifs","%","% cumules")
  res = res[,c(1,indic+1)]
  res=as.data.frame(res)
  
  return(res)
  
}

### Gerate 1D plots
generateUniPlot<-function(input,data)
{
  
  ind = which(colnames(data)%in%input$UniVar)
  
  dataTmp = data[,ind]
  TestNum = is.numeric(dataTmp)
  dataTmp = data.frame(x=dataTmp)
  namesTmp = names(data)[ind]
  gg = NULL
  
  ## Numeric data
  if(TestNum && !is.null(input$RadioPlotUni))
  {
    if(input$RadioPlotUni=="hist")
    {
      gg = ggplot(dataTmp,aes(x=x))  + xlab(namesTmp) + theme_bw()
      if(input$HistDens=="freq") gg = gg + geom_histogram(binwidth = input$binwidth,size=input$SizeQQplot-0.5,fill=input$ColorUniplot,alpha=input$TransAlphaUni/100, color="black",aes(y = ..density..))
      if(input$HistDens=="counts") gg = gg + geom_histogram(binwidth = input$binwidth,size=input$SizeQQplot-0.5,fill=input$ColorUniplot,alpha=input$TransAlphaUni/100, color="black")
      if(input$CheckDens) gg = gg + geom_density(size=input$SizeQQplot-0.5)
      if(input$SensGraph=="Hori") gg = gg + coord_flip()  
    }
    
    if(input$RadioPlotUni=="box")
    {
      gg = ggplot(dataTmp,aes(1,x))  + geom_boxplot(fill=input$ColorUniplot,alpha=input$TransAlphaUni/100,size=input$SizeQQplot)+xlim(c(0,2)) 
      gg = gg + ylab(namesTmp) + xlab("")+ theme_bw()
      if(input$CheckAddPointsBox) gg = gg + geom_jitter()
      if(input$SensGraph=="Hori") gg = gg + coord_flip()
      
    }
    
    if(input$RadioPlotUni=="densities")
    {
      gg = ggplot(dataTmp,aes(x=x))  + theme_bw()
      gg = gg + geom_density(fill=input$ColorUniplot,alpha=input$TransAlphaUni/100,size=input$SizeQQplot-0.5) + xlab(namesTmp)
      if(input$SensGraph=="Hori") gg = gg + coord_flip() 
      
    }
    
    if(input$RadioPlotUni=="qqplot")
    {
      gg = ggplot(dataTmp, aes(sample = x)) + geom_point(stat = "qq",color = input$ColorUniplot,alpha=input$TransAlphaUni/100,size=input$SizeQQplot)
      gg = gg + ylab(namesTmp)+ theme_bw()
    }
    
    return(gg)
  }
  
  # Quali data
  if(!TestNum && !is.null(input$RadioPlotUni))
  {
    if(input$RadioPlotUni=="BarPlot")
    {
      gg = ggplot(dataTmp, aes(x))  + xlab(namesTmp)+ theme_bw()
      gg = gg+  geom_bar(fill=input$ColorUniplot,alpha=input$TransAlphaUni/100,size=input$SizeQQplot,  color="black",width=input$widthBarPlot/100)
      if(input$SensGraph=="Hori") gg = gg + coord_flip()
      if(input$BarCircular)  gg = gg + coord_polar()
    }
    
    if(input$RadioPlotUni=="Pie")
    {
      
      count = table(dataTmp$x)
      dataTmp2 = data.frame(frac = count/sum(count), xUnique = as.factor(names(count)))
      dataTmp2 = dataTmp2[order(dataTmp2$frac.Freq), ]
      dataTmp2$ymax = cumsum(dataTmp2$frac.Freq)
      dataTmp2$ymin = c(0, head(dataTmp2$ymax, n=-1))
      
      dataTmp2$xminPie = 1-input$PieWidth/100
      gg =  ggplot(dataTmp2, aes(fill=xUnique, ymax=ymax, ymin=ymin, xmax=4, xmin=4*xminPie)) + geom_rect(alpha=input$TransAlphaUni/100) + coord_polar(theta="y") 
      gg =  gg + xlim(c(0, 4)) + scale_fill_discrete(name=namesTmp)+ theme_bw()
      
    }
    return(gg)
  }
  
  
  
}


### Gerate 2D plots
generateBiPlot<-function(input,rangesBiplot,data)
{
  
  var1 = input$VariableSelectBi1
  var2 = input$VariableSelectBi2
  
  ind1  = which(colnames(data)%in%c(var1))
  ind2  = which(colnames(data)%in%c(var2))
  ind =  c(ind1,ind2)
  
  Num = sapply(data[,unique(ind)],is.numeric)
  
  if(length(unique(ind))==2)
  {
    if(length(which(Num))==2) data2 = data.frame(x=data[,ind1],y=data[,ind2])
    if(length(which(Num))==1) data2 = data.frame(x=data[,ind[which(Num)]],y=data[,ind[which(!Num)]])        
  }
  
  if(length(unique(ind))==1)  data2 = data.frame(x=data[,ind],y=data[,ind])
  
  if(!is.null(input$RadioPlotBi))
  {
    if(input$RadioPlotBi=="Nuage") 
    {
      gg = ggplot(data2,aes(x,y))  
      gg = gg + geom_point(color=input$ColorBiplot,size=input$SizePoint,alpha=input$TransAlphaBi/100) + theme_bw()
      if(!is.null(rangesBiplot$x) && !is.null(rangesBiplot$x)) gg = gg + xlim(rangesBiplot$x) + ylim(rangesBiplot$y)
      gg = gg + xlab(var1) + ylab(var2)
      if(input$CheckLM) gg = gg + geom_smooth(method='lm')
    }    
    
    if(input$RadioPlotBi=="densities" )
    {
      if(length(which(Num))==2)
      {
        vartmp = c(data[,ind[1]],data[,ind[2]])
        y=c(rep(var1,nrow(data)),rep(var2,nrow(data)))
        data2 = data.frame(x=vartmp,y=y)
      }
      
      gg = ggplot(data2,aes(x=x,fill=y)) 
      gg = gg + xlab("Valeurs")+labs(colour = var2) + scale_fill_discrete(name="Légende")
      gg = gg + geom_density(size=input$SizePoint,alpha = input$TransAlphaBi/100) + theme_bw()
      if(input$SensGraphBi=="Hori") gg = gg + coord_flip()
    }    
    
    if(input$RadioPlotBi=="hist" )
    {
      if(length(which(Num))==2)
      {
        vartmp = c(data[,ind[1]],data[,ind[2]])
        y=c(rep(var1,nrow(data)),rep(var2,nrow(data)))
        data2 = data.frame(x=vartmp,y=y)
      }
      
      gg = ggplot(data2,aes(x=x,fill=y)) 
      gg = gg + xlab("Valeurs")+labs(colour = var2) + scale_fill_discrete(name="Légende")
      gg = gg + geom_histogram(binwidth = input$binwidthBi,size=input$SizePoint-0.5,alpha = input$TransAlphaBi/100,color="black") + theme_bw()
      if(input$SensGraphBi=="Hori") gg = gg + coord_flip()
    }    
    
    if(input$RadioPlotBi=="box" )
    {
      if(length(which(Num))==2)
      {
        vartmp = c(data[,ind[1]],data[,ind[2]])
        y=c(rep(var1,nrow(data)),rep(var2,nrow(data)))
        data2 = data.frame(x=vartmp,y=y)
      }
      labtmp = unique(data2$y)
      gg = ggplot(data2,aes(y,x)) + scale_fill_discrete(name="Légende")
      gg = gg+ geom_boxplot(aes(fill=y),alpha=input$TransAlphaBi/100,size=input$SizePoint-0.5) + xlab("")+ theme_bw()
      if(input$CheckAddPointsBoxBi) gg = gg + geom_jitter(size=input$SizePoint)
      if(input$SensGraphBi=="Hori") gg = gg + coord_flip()
      
    }    
    
    return(gg)
    
  }
  
}


