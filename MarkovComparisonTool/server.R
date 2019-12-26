#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
require(expm)
require(pracma)
require(tidyverse)
require(matrixcalc)



require(expm)
require(pracma)
require(tidyverse)


#The engine for the Markov model conversion
matroot<-function(A,p)
{
  if(is.integer(1/p))
  {
    return(A%^%(1/p))
  }
  else if (!is.singular.matrix(eigen(A)$vectors))
  {
    return(eigen(A)$vectors %*% diag(eigen(A)$values^{1/p}) %*% solve(eigen(A)$vectors))
  }
  else {return(matrix(0,1,1))}
  }


is.stochastic<-function(A=diag(1,nrow=2,ncol=2))
{
  rowsums<-round(rowSums(A),15)
  rowmins<-numeric()
  if(any(Im(A)!=0))
  {return(FALSE)}
  else{
    
    if (all(rowsums==1) & min(A) >= 0 )
    {return(TRUE)}
    else {return(FALSE)}
  }}

markovTest<-function(oldM = matrix(c(0.7,0.3,0,0,0.2,0.8,0,0,1),nrow=3,byrow=T),
                     newM = matrix(c(0.97,0.03,0,0,0.87,0.13,0,0,1),nrow=3,byrow=T)
                     ,oldCycle=12,newCycle=1,initDist=0,ncycles=10){
  
  Mrow=nrow(oldM)
  Mcol=ncol(oldM)
  # If one matrix is not square or if one is a different size than the other then return an error and exit.
  
  #if (initDist==0)
  #{
  #  initDist<-c(1,rep(0,Mcol-1))
  #}
  
  
  # Cycle length and transition matrices to use for output
  
  if ( is.stochastic(matroot(newM,newCycle/oldCycle)))
  {
    outCycle<-oldCycle
    baseM<-oldM
    altM<-matroot(newM,newCycle/oldCycle)
    print(paste0("Using old cycle length ",oldCycle))
  }
  else if ( is.stochastic(matroot(oldM,oldCycle/newCycle)))
  {
    outCycle<-newCycle
    altM<-newM
    baseM<-matroot(oldM,oldCycle/newCycle)
    print(paste0("Using new cycle length ",newCycle))
  }
  else
  {
    outCycle<-Lcm(oldCycle,newCycle)
    baseM<-matroot(oldM,oldCycle/outCycle)
    altM<-matroot(newM,newCycle/outCycle)
    print(paste0("Using LCM of cycle lengths ",outCycle))
  }
  
  
  
  
  #cohort=data.frame(initDist)
  oldDist<-matrix(initDist,nrow=1)
  newDist<-matrix(initDist,nrow=1)
  errors<-numeric()
  avg<-data.frame()
  
  for ( i in 1:ncycles)
  {
    oldDist<-rbind(oldDist,oldDist[i,]%*%(baseM) )
    newDist<-rbind(newDist,newDist[i,]%*%(altM) )
    errors<-rbind(errors,sum(abs(oldDist[i+1,]-newDist[i+1,]))/2)
  }
  
  #print("Cohort distribution per cycle using original transition matrix:")
  #print(oldDist)
  #print("Cohort distribution using new transition matrix:")
  #print(newDist)
  #print("Proportion of Cohort Misclassified by cycle:")
  #print(errors)
  
  outDist=list(oldDist=data.frame(cbind((0:ncycles)*outCycle),oldDist),
               newDist=data.frame(cbind((0:ncycles)*outCycle),newDist))
  names(outDist$oldDist)<-c("TimeElapsed",paste0("State",1:(ncol(outDist$oldDist)-1)))
  names(outDist$newDist)<-c("TimeElapsed",paste0("State",1:(ncol(outDist$newDist)-1)))
  outDist$oldDist$sourceM<-"Original"
  outDist$newDist$sourceM<-"Converted"
  return(outDist)
}





# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$nstates<-reactive({input$nstates})
  FoldM<-reactive({as.matrix(hot_to_r(input$oldM))})
  FnewM<-reactive({as.matrix(hot_to_r(input$newM)) })
  FconvM<-reactive(matroot(FoldM(),(input$oldCycle/input$newCycle)))
  FconvM2<-reactive(matroot(FnewM(),(input$newCycle/input$oldCycle)))
  FinitD<-reactive({as.matrix(hot_to_r(input$initD))})
 
   
   output$summary1 <- renderPrint({
     
     if (is.stochastic(FconvM()))
     {
       return(paste0("In this case, the transition matrix A with cycle length ",
                     input$oldCycle,
              "can be converted to a transition matrix with cycle length ",
              input$newCycle,
              ". The corresponding matrix is:"))
            }
     else{
       return(
         paste0("In this case, the transition matrix A with cycle length ",
                input$oldCycle,
                "cannot be converted to a transition matrix with cycle length ",
                input$newCycle,
                ". The attempt to convert A results in a non-stochastic matrix (i.e. a matrix with negative and/or complex entries, which cannot be probabilities):")
         
       )
     }
   })

   
   
    output$convM<-renderTable({
                    FconvM()
                               
      })
     
    output$summary2<-renderPrint(
      {
        if ( is.stochastic(FconvM()) & all(round(FconvM(),16)==round(FnewM(),16)))
          {return("This is the same as the matrix B, so there is no error introduced in converting the cycle length.")}
        else if ( is.stochastic(FconvM()) & !all(round(FconvM(),16)==round(FnewM(),16)))
          {return("This is not the same as the matrix B.  The matrix shown above is the correct conversion of A and should be used instead. The plot below shows the proportion of the cohort that are incorrectly classified when using the matrix B.")}
          
          {
          return("As a consequence, the new matrix B does not represent the same transition probabilities as A over time. The plot below shows the proportion of the cohort that are incorrectly classified when using the matrix B.")
        }
        
        }
      
    )
    
    output$summary3<-renderText({
      if (is.stochastic(FconvM2()))
      {
        return(paste0("In this case, the transition matrix B with cycle length ",
                      input$newCycle,
                      " can be converted to a transition matrix with cycle length ",
                      input$oldCycle,
                      ", so the two matrices will be compared using the latter cycle length."))
      }
      else if (!is.stochastic(FconvM2()) & is.stochastic(FconvM()))
      {
        return(paste0("In this case, the transition matrix B with cycle length ",
                      input$newCycle,
                      " cannot be converted to a transition matrix with cycle length ",
                      input$oldCycle,
                      ". However, the matrix A can be converted to a matrix with cycle length",
                      input$newCycle,
                      ", so the two matrices will be conpared using this cycle length."))
      }      
      if (!is.stochastic(FconvM2()) & is.stochastic(FconvM()))
      {
        return(paste0("In this case, the cycle length of A cannot be converted to the cycle length of B, nor can the cycle length of B be converted to that of A.  For this reason, the matrices will be compared using the smallest possible common cycle length of",
               Lcm(input$oldCycle,input$newCycle)
               ))
      }
   } )
    
   
  output$summary4<-renderText({
    if (is.stochastic(FconvM2()) & !is.stochastic(FconvM()))
    {
      sprintf("Note: it is possible to convert B exactly to a matrix with cycle length %i.  Is this an option in your model? ",input$oldCycle)
    }
    
  })  
    
   output$oldM <- renderRHandsontable({
     defaultoldM<- 
       matrix(
         if ( input$nstates ==3 )
         {
           matrix(c(0.7,0.3,0,0,0.2,0.8,0,0,1),nrow=3,byrow=T)
         }
         else
         {diag(1,input$nstates)}
         
         , nrow=input$nstates
       )
    
    
    
     rhandsontable(data.frame(defaultoldM),readOnly = FALSE)
   })
   
   output$errorPlot<-renderPlot({
     errordf<- select(MarkovDistDF(),
                      Time.Elapsed,
                      Proportion.Misclassified
     )  
     return(ggplot(data=errordf,aes(x=Time.Elapsed,y=Proportion.Misclassified))+geom_line(size=1.5)+labs(title="Proportion of cohort misclassified over time"))
     
   })
   
   
   
   output$cohortPlot<-renderPlot({
     mdf<-rbind(gather(MTest()$oldDist,key=state,value=proportion,-sourceM,-TimeElapsed),
                gather(MTest()$newDist,key=state,value=proportion,-sourceM,-TimeElapsed))
     
     ggplot(data=mdf,aes(x=TimeElapsed,y=proportion,colour=state,fill=state))+
       geom_col(position = "stack")+facet_wrap(sourceM~.,ncol = 1)
     
   })
   
   output$initD<-renderRHandsontable(
     {
       rhandsontable(data.frame(t(c(1,rep(0,input$nstates-1)))))
     }
     
   )
   output$newM<-renderRHandsontable({
     defaultnewM<- 
       matrix(
         if ( input$nstates ==3 )
         {
           matrix(c(0.97,0.03,0,0,0.87,0.13,0,0,1),nrow=3,byrow=T)
         }
         else
         {diag(1,input$nstates)}
         
         , nrow=input$nstates
       )
     
     
     rhandsontable(data.frame(defaultnewM),readOnly = FALSE)
   })
   
   
   MTest<-reactive({
            markovTest(FoldM(),
                       FnewM(),
                       input$oldCycle,
                       input$newCycle,
                       FinitD(),
                       input$ncycles)
       })
   
   
   MarkovDistDF<-reactive(
      {
       mtest<-MTest()
       diffs<-rowSums(abs(select(mtest$oldDist,-sourceM)-select(mtest$newDist,-sourceM)))/2
       names(mtest$oldDist)<-paste0("Old.",names(mtest$oldDist))
       names(mtest$newDist)<-paste0("New.",names(mtest$newDist))
       mtest2<-data.frame(mtest$oldDist,mtest$newDist,"Proportion.Misclassified"=diffs)%>%                                                    
         select(-New.TimeElapsed,
                -Old.sourceM,
                -New.sourceM)
       names(mtest2)[1]<-"Time.Elapsed"
       return(mtest2)
       
     }
     )
   
   output$MarkovDist<-renderRHandsontable({
                                rhandsontable(MarkovDistDF(),readOnly = TRUE)    
   })
   
   
   
   myHeightAlgorithm<-reactive({
   
   400*input$nstates
   
   })
   
   output$plots<-renderPlot({
     
       mdf<-rbind(gather(MTest()$oldDist,key=state,value=proportion,-sourceM,-TimeElapsed),
                gather(MTest()$newDist,key=state,value=proportion,-sourceM,-TimeElapsed))
     
       ggplot(data=mdf,aes(x=TimeElapsed,y=proportion,colour=sourceM,fill=sourceM))+
       geom_col(position = "dodge")+facet_wrap(state~.,ncol = 1)
       #theme(aspect.ratio = 1/input$nstates)
     
      }, height=myHeightAlgorithm)
   
   })
