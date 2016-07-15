
#1) For now, sample name in reference file must be exactly "Father", "Mother", "Child" if one uses built in pedigrees.
#2) Should we have more options for pedigrees? For now, if custom pedigree is specified for one
#pedigree, custom must also be specified for the other pedigree
#3) Files must have '.' as decimal separator
#4) Change dropout windows to gedit?
#5) Contributors must be the same under both hypotheses. Their familiy relationship can be changed with
#   the pedigree, but any reference profiles must be the same under both hypotheses

relMixGUI <- function(){
  
  options("guiToolkit"="tcltk")

  #Files must have '.' as decimal separator
  tableReader <- function(filename) { #helpfunction
    tab <- read.table(filename,header=TRUE,sep="\t",stringsAsFactors=FALSE,na.strings=c(NA,""))
    tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings=c(NA,"")) } ,error=function(e) e)
    tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE,na.strings=c(NA,"")) } ,error=function(e) e)
    if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE,na.strings=c(NA,""))
    return(tab) #need dataframe to keep allele-names correct!!
  }
  
  tableWriter <- function(filename,obj){
    write.table(obj,file=filename,sep="\t",row.names=FALSE,quote=FALSE)
  }

  f_importprof <- function(h,...) {
    type=h$action #get type of profile
    proffile = gfile(text=paste("Open ",type," file",sep=""),type="open",
                     filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*"))))
    Data <- tableReader(proffile) #load profile
    assign(h$action,Data,envir=mmTK) #save object
  }
  
  f_export <- function(obj) {
    savefile = gfile(text=paste("Save file as",sep=""),type="save", initialfilename = "LR.txt")
                     #filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*")))
    tableWriter(savefile,obj) #load profile
  }

  #Make pedigree
  f_pedigree <- function(h,...){
    if(svalue(h$obj)=="Paternity"){
      #Define the persons involved in the case
      persons <- c("Mother", "Father", "Child")
      sex <- c("female", "male", "male")
      ped1 <- FamiliasPedigree(id=persons, dadid=c(NA,NA,"Father"), momid=c(NA,NA,"Mother"), sex=c("female", "male", "male"))
    } else if(svalue(h$obj)=="Unrelated"){
      #Define the persons involved in the case
      persons <- c("Mother", "Father", "Child")
      sex <- c("female", "male", "male")
      ped1 <- FamiliasPedigree(id=persons, dadid=c(NA,NA,NA), momid=c(NA,NA,"Mother"), sex=c("female", "male", "male"))
    } else {
      pedfile = gfile(text=paste("Open pedigree R file",sep=""),type="open",
                      filter=list("text"=list(patterns=list("*.R")),"all"=list(patterns=list("*"))))
      source(pedfile)
      if(!all(exists("persons"),exists("ped1"))) stop("File should define both pedigree and persons")
    }
    assign(paste("ped",h$action,sep=""),ped1,envir=mmTK)
    assign(paste("persons_ped",h$action,sep=""),persons,envir=mmTK)
  }
  
  # f_pedigree <- function(h,...){
  #   
  #   if(svalue(h$obj)=="Paternity"){
  #     #Define the persons involved in the case
  #     #If pedigree 2, use same individuals as in first pedigree
  #     if(svalue(h$action)=='2'){
  #       firstPed <- get('ped1',envir=mmTK)
  #       persons1 <- get('persons_ped1',envir=mmTK)
  #       persons <- persons1
  #     } else{
  #       persons <- c("Mother", "Father", "Child")
  #       sex <- c("female", "male", "male")
  #     }
  #     #persons <- c("Mother", "Father", "Child")
  #     #sex <- c("female", "male", "male")
  #     ped1 <- FamiliasPedigree(id=persons, dadid=c(NA,NA,"Father"), momid=c(NA,NA,"Mother"), sex=c("female", "male", "male"))
  #   } else if(svalue(h$obj)=="Unrelated"){
  #     #Define the persons involved in the case
  #     if(svalue(h$action)=='2'){
  #       firstPed <- get('ped1',envir=mmTK)
  #       persons1 <- get('persons_ped1',envir=mmTK)
  #       persons <- persons1
  #     } else{
  #       persons <- c("Mother", "Father", "Child")
  #       sex <- c("female", "male", "male")
  #     }
  #     #persons <- c("Mother", "Father", "Child")
  #     #sex <- c("female", "male", "male")
  #     ped1 <- FamiliasPedigree(id=persons, dadid=c(NA,NA,NA), momid=c(NA,NA,"Mother"), sex=c("female", "male", "male"))
  #   } else {
  #     pedfile = gfile(text=paste("Open pedigree R file",sep=""),type="open",
  #                     filter=list("text"=list(patterns=list("*.R")),"all"=list(patterns=list("*"))))
  #     source(pedfile)
  #     if(!all(exists("persons"),exists("ped1"))) f_errorWindow("File should define both pedigree and persons")#stop("File should define both pedigree and persons")
  #   }
  #   assign(paste("ped",h$action,sep=""),ped1,envir=mmTK)
  #   assign(paste("persons_ped",h$action,sep=""),persons,envir=mmTK)
  # }

  #Get values from object and assign to new object in mmTK environment
  f_values <- function(h,...) {
    r <- svalue(h$obj)
    assign(h$action,r,envir=mmTK)
  }
  #Pop-up error window
  f_errorWindow <- function(message){
    errorWindow <- gwindow("Error")
    glabel(message,container=errorWindow)
    gbutton("ok", container=errorWindow,handler=function(h,...){
      dispose(h$obj)
    })
  }

  #Makes pop-up window for database
  f_database <- function(h,...){
    freqWindow <- gwindow("Frequency database",visible=FALSE)
    freqGroup2 <- ggroup(horizontal = TRUE, container=freqWindow,spacing=7)
    gbutton(text="Import allele frequencies",container=freqGroup2,handler=f_importprof,action="frequencies")
    optButton <- gbutton("Options", handler=f_options, container=freqGroup2)
    saveButton <- gbutton(text = 'ok',container=freqWindow, handler = function(h,...){
      
      #Get frequencies
      if(!exists('frequencies',envir=mmTK)) {f_errorWindow("Alle frequencies not imported")}
      if(!exists('mixture',envir=mmTK)) {f_errorWindow("Import mixture profile before frequencies")}
      if(!exists('reference',envir=mmTK)) {f_errorWindow("Import reference profiles before frequencies")}
      A <- get("frequencies",envir=mmTK)
      M <- get("mixture",envir=mmTK)
      G <- get("reference",envir=mmTK)
      
      freqs <- A[,-1,drop=FALSE]
      colnames(freqs) <- colnames(A)[2:ncol(A)]
      rownames(freqs) <- A[,1]
      optPar <- get('optPar',envir=mmTK)
      
      freqsS <- f_silent(freqs,optPar$silent) #Add silent allele
      freqsS <- get('freqsS',envir=mmTK)
      db2 <- f_unobserved(freqsS,M,G,optPar$MAF) #Add unobserved alleles
      
      dispose(freqWindow)
    }) #end handler savebutton
    
    visible(freqWindow) <- TRUE
  }


  #Makes pop-up window to fill in mutation details
  f_mutations <- function(h,...){
    mutWindow <- gwindow("Mutations",visible=FALSE)
    mutGroup2 <- ggroup(horizontal = TRUE, container=mutWindow,spacing=7)
    mutFrame2 <- gframe("Mutation model",container=mutGroup2)
    mutFrame3 <- gframe("Range",container=mutGroup2)
    #Range window
    objRange <- gedit('0.5', width=4, container=mutFrame3)
    enabled(objRange) <- FALSE
    #Common male and female mutation model
    objMutModel <- gcombobox(c("Equal","Proportional","Stepwise"), selected=1,
                              container=mutFrame2, handler=function(h,...) {
                               assign("Mut",svalue(h$obj),envir=mmTK)
                               #Range will not be used unless model is Stepwise
                               assign("range1",0.5,envir=mmTK)
                               f_changeMutation(svalue(h$obj),objRange)})
    #List to store parameters in
    objMut <- vector('list',3)
    names(objMut) <- c("mutRange","fMutRate","mMutRate")
    objMut[[1]] <- objRange
    mutFrame4 <- gframe("Mutation rates",container=mutWindow,horizontal=TRUE)
    #Separate male and female mutation rates
    glabel("Female",container=mutFrame4)
    objMut[[2]] <- gedit('0', width=5, container=mutFrame4)
    glabel("Male",container=mutFrame4)
    objMut[[3]] <- gedit('0', width=5,container=mutFrame4)
    saveButton <- gbutton(text = "Save", handler = function(h,...) saveParameters(objMut,mutWindow,"mutPar"),container=mutWindow)
    visible(mutWindow) <- TRUE
  }

  #Set mutation model
  f_changeMutation <- function(mutMod,objRange){
    if(mutMod=="Stepwise"){ #If stepwise, also give option to give a mutation range
      #Common male and female mutation range
      enabled(objRange) <- TRUE
    }
    if(mutMod=="Equal"){
      enabled(objRange) <- FALSE
    }
    if(mutMod=="Proportional"){
      enabled(objRange) <- FALSE
    }
    #if(mutMod=="Custom"){ #If custom, read in mutation model from file
    #  enabled(objRange) <- FALSE
    #  enable(objRateF)<- FALSE
    #  enable(objRateM) <- FALSE
    #  mutfile = gfile("Open mutation R file",type="open",
    #                  filter=list("text"=list(patterns=list("*.R",".txt")),"all"=list(patterns=list("*"))))
    #  MM <- read.table(mutfile)
    #  assign(paste(h$action,"Mat",sep=""),MM,envir=mmTK)
    #}
  }

  #Makes option pop-up window to fill in theta and silent alleles
  f_options <- function(h,...){
    optWindow <- gwindow("Options",visible=FALSE)
    optGroup <- ggroup(horizontal = TRUE, container=optWindow,spacing=7)
    #Theta
    glabel("Theta",container=optGroup)
    objOpt <- vector('list',3)
    names(objOpt) <- c("theta","silent","MAF")
    objOpt[[1]] <- gedit('0', width=5, container=optGroup)
    #Silent allele
    glabel("Silent allele frequency",container=optGroup)
    objOpt[[2]] <- gedit('0', width=5, container=optGroup)
    #Minimum allele frequency
    glabel("Min. allele frequency",container=optGroup)
    objOpt[[3]] <- gedit('0.001', width=5, container=optGroup)
    saveButton <- gbutton(text = "ok", ,container=optWindow, handler = function(h,...){
      saveParameters(objOpt,optWindow,"optPar")
    })
    visible(optWindow) <- TRUE
  }

  #Makes pop-up window to fill in dropout values per contributor and drop-in
  f_dropout <- function(h,...){
    dropWindow <- gwindow("Dropout and drop-in",visible=FALSE)
    #Get pedigree data
    if(!exists("idxC",envir=mmTK)){
      f_errorWindow("Specify contributors first")
    }
    dropFrame2 <- gframe("Dropout for each contributor",container=dropWindow)
    #Get contributors
    idxC <- get("idxC",envir=mmTK)
    objDrop <- list()
    for(i in 1:length(idxC)){
      g <- glabel(idxC[i], container=dropFrame2,horizontal=FALSE)
      objDrop[[i]] <- gspinbutton(0,1,by=0.01, value=0,container=dropFrame2)
    }
    dropinFrame <- gframe("Drop-in",container=dropWindow)
    objDrop[[length(idxC)+1]] <- gspinbutton(0,1,by=0.01, value=0,digits=3,container=dropinFrame)
    names(objDrop) <- c(idxC,"dropin")
    saveButton <- gbutton(text = "ok",container=dropWindow, handler = function(h,...) {
      saveParameters(objDrop,dropWindow,"dropliste")
      })
    visible(dropWindow) <- TRUE
  }
  #Save button
  saveParameters <- function(objList,window,varName){
    values <- list()
    for(i in 1:length(objList)){
      #If comma used as decimal, replace with dot
      v <- as.numeric(sub(",",".",svalue(objList[[i]])))
      if(length(v)==0) {f_errorWindow("Specify all values"); stop()}
      else if(is.na(v) | (v<0 | v>1) ) {f_errorWindow("Need values in the range 0.0 and 1.0"); stop()}
      else values[[i]] <- v
    }
    names(values) <- names(objList)
    assign(varName,values,envir=mmTK)
    dispose(window)
  }

  f_contributors <- function(h,...){
    #Get pedigree data
    if(!exists("persons_ped1",envir=mmTK)){
      f_errorWindow("Specify pedigrees first")
    }
    persons1 <- get("persons_ped1",envir=mmTK)
    persons2 <- get("persons_ped2",envir=mmTK)
    if(!identical(persons1,persons2)) {
      f_errorWindow("Persons in pedigree 1 and 2 must match!")
    } else{
      contWindow <- gwindow("Contributors")
      contFrame <- gframe("Specify contributors in mixture",container=contWindow)
      gcheckboxgroup(union(persons1,persons2), checked=FALSE,horizontal=FALSE,container=contFrame,
                     handler=f_values,action="idxC")
   
      
      # #Different contributors under each hypothesis
      # contFrame <- gframe("Contributors in pedigree 1",container=contWindow)
      # gcheckboxgroup(union(persons1,persons2), checked=FALSE,horizontal=FALSE,container=contFrame,
      #                handler=f_values,action="idxC1")
      # contFrame2 <- gframe("Contributors in pedigree 2",container=contWindow)
      # gcheckboxgroup(union(persons1,persons2), checked=FALSE,horizontal=FALSE,container=contFrame2,
      #                handler=f_values,action="idxC2")
      
      gbutton("ok", container=contWindow,handler=function(h,...){
        dispose(h$obj)
      })
    }
   }

  #Get mixture in the right format
  f_mixture <- function(E){
    m <- length(unique(E$Marker)) #Number of markers
    #n <- length(unique(E$SampleName)) #Number of samples (if replicates)
    mix <- split(E[,-c(1,2)],E$Marker) #Split according to markers
    #Remove NA's
    lapply(mix,function(x) x[!is.na(x)])
    #lapply(mix,function(x) { #Split according to sample (if replicates)
    #  x <- x[,!is.na(x),drop=F]
    #  split(x[-c(1,2)],x$SampleName)
    #})
  }
  #Reads genotypes of known contributors and stores in list
  f_genotypes <- function(GT){
    m <- length(unique(GT$Marker)) #Number of markers
    n <- length(unique(GT$SampleName)) #Number of contributors
    gt <- split(GT,GT$Marker) #Split according to markers
    lapply(gt,function(x) split(x[,3:4],x$SampleName)) #Split according to individual
  }
 
  #Add silent allele to database
  f_silent <- function(freqs,ps){
    
    if(ps>0) {
      freqsS <- rbind(freqs,rep(as.numeric(ps),ncol(freqs)))
      rownames(freqsS) <- c(rownames(freqs),'Silent')
    } else freqsS <- freqs
    assign('freqsS',freqsS,envir=mmTK)
  }
  
  #Add alleles not in database with frequency MAF
  #Checks if any alleles are below MAF (by running f_MAF)
  f_unobserved <- function(freqs,M,G,MAF){
    
    alleleNames <- rownames(freqs)
    markerNames <- colnames(freqs)
    
    n <- ncol(M)
    keepIx <- alNotDB <- numeric()
    #Go through each marker in mixture profile to look for alleles not in database
    #Simultaneously checking in reference profile (assumes that markers in reference and
    #mixture are the same)
    for(i in 1:nrow(M)){
      #Get all database frequencies for the marker
      mark <- M[i,2]
      ix <- which(mark==markerNames)
      al <- alleleNames[!is.na(freqs[,ix])]
      #Alleles in mixture
      am <- M[i,3:n][!is.na(M[i,3:n])] 
      #Alleles in genotypes
      ag <- unlist(G[G[,2]==mark,3:4])
      aa <- unique(c(am,ag))
      #Check if all alleles in mixture and reference exist in database
      idx <- c(aa%in%al)
      
      if(any(!idx)){
        keepIx <- c(keepIx,rep(ix,sum(!idx))) #Index of marker
        alNotDB <- c(alNotDB,aa[!idx]) #Alleles not found in db
      }
    }
    
    #Change format of database before adding new alleles
    db <- numeric()
    for(i in 1:ncol(freqs)){
      ix <- which(!is.na(freqs[,i]))
      a <- alleleNames[ix]
      f <- freqs[ix,i]
      Marker <- rep(markerNames[i],length(ix))
      dbNew <- data.frame(Marker,a,f)
      colnames(dbNew) <- c("Marker","Allele","Frequency")
      db <- rbind(db,dbNew)
    }
    
    if(length(alNotDB)>0) { #There are alleles not in database
      
      mess <- paste("Allele",alNotDB, "will be added to marker",markerNames[keepIx], "with frequency",MAF)
      gmessage(mess, title="Note",icon = "info")
      # gconfirm(mess, title="Note",icon = "info",handler=function(h,...){ #Begin handler add alleles
         #Add new allele at the end of database
         newData <- data.frame(markerNames[keepIx],alNotDB,MAF)
         colnames(newData) <- c('Marker','Allele','Frequency')
         db <- data.frame(Marker=c(as.character(db$Marker),as.character(newData$Marker)),
                          Allele=c(as.character(db$Allele),as.character(newData$Allele)),
                          Frequency=c(db$Frequency,newData$Frequency))
         assign('db',db,mmTK)
         #Check for MAF, scale and sort
         f_MAF(db,MAF)
      # }) #end handler add alleles
      
    } else{ #No alleles added to database
      #Check for MAF, scale and sort
      assign('db',db,mmTK)
      f_MAF(db,MAF)
    }
  }# end f_unobserved
  
  #Check if any allele frequencies are below MAF and sets frequency to MAF
  #Scales frequencies if necessary (with f_MAF)
  #Function used by f_unobserved
  f_MAF <- function(db,MAF){
    db <- get('db',mmTK)
    if(any(db$Frequency<MAF)){ #Frequencies below MAF
      gconfirm("Some frequencies are below the min. allele frequency.
          Change the indicated frequencies?", title="Note",icon = "question",handler = function(h,...){
        db$Frequency[db$Frequency<MAF] <- MAF
        assign('db',db,mmTK)
      })
      
    }
    #Check if scaling is necessary
    db <- get('db',mmTK)
    f_scale(db)
  }
  
  
  #Check that all frequencies sum to 1, otherwise scale
  #Sort database and assign final database to environment
  #Function used by f_MAF
  f_scale <- function(db){
    db <- get('db',mmTK)
    markerNames <- unique(db[,1])
    
    #Sort database according to marker, then allele. Silent allele last
    #First reorder levels of Allele
    aL <- levels(db[,2])
    if(any(aL=='Silent')) { db$Allele <- factor(db$Allele,c(sort(as.numeric(aL[which(!aL=='Silent')])),'Silent'))
    } else db$Allele <- factor(db$Allele,sort(as.numeric(aL)))
    
    db <- db[order(db$Marker,db$Allele),]
    #Assign database that will be used if scaling is not done
    assign('dbF',db,envir=mmTK) #Final database
    
    #Check that frequencies sum to 1, otherwise scale
    #Round to 4 decimals
    sums <- sapply(1:length(markerNames),function(i) round(sum(db[db[,1]==markerNames[i],3]),4))
    ix <- which(sums!=1)
    if(length(ix)>0) { 
      #TO DO: Add rest allele if no scaling. What do we do if sum > 1?
      #Add rest allele
      #mSums <- markerNames[ix]
      #dbF <- rbind(db,data.frame(Marker=mSums,Allele='r',Frequency=1-sums[ix]))
      #assign('dbF',dbF,envir=mmTK) #Final database if no scaling
      gconfirm("Frequencies do not sum to 1. Do you want to scale?", title="Note",icon = "info",handler = function(h,...){
        mSums <- markerNames[ix]
        for(m in mSums){
          db[db[,1]==m,3] <- db[db[,1]==m,3]/sum(db[db[,1]==m,3])
        }
        assign('dbF',db,envir=mmTK) #Final database
      })
    }
  }
  
    

  f_LR <- function(){
    
    ####### Get input ########
    #Mutation parameters
    mutPar <- get("mutPar",envir=mmTK)
    r1 <- mutPar$fMutRate
    r2 <- mutPar$mMutRate
    range1 <- mutPar$mutRange
    mutModel <- get("Mut",envir=mmTK)
    optPar <- get('optPar',mmTK)
    theta <- optPar$theta

    #Profiles
    if(!exists('mixture',envir=mmTK)){ f_errorWindow("Mixture not imported"); stop()}
    if(!exists('reference',envir=mmTK)){ f_errorWindow("Reference profiles not imported"); stop()}
    if(!exists('frequencies',envir=mmTK)){ f_errorWindow("Missing allele frequencies"); stop()}
    if(!exists("ped1",envir=mmTK)){ f_errorWindow("Missing pedigree information"); stop()}
    if(!exists("ped2",envir=mmTK)){ f_errorWindow("Missing pedigree information"); stop()}
    if(!exists("idxC",envir=mmTK)){ f_errorWindow("Contributors not specified"); stop()}
    
    E <- get("mixture",envir=mmTK) #get object
    G <- get("reference",envir=mmTK) #get object
    #Remove AMEL marker? Or not allow for it?
    R <- f_mixture(E)
    knownGenos <- f_genotypes(G)
    #Frequencies
    db2 <- get("dbF",envir=mmTK) 
    alleleNames <- as.character(db2[,2])
    
    #Check if there are non-numeric allele names and stepwise mutation model
    #If so, give a warning
    isNumeric <- suppressWarnings(as.numeric(alleleNames[-which(alleleNames=="Silent")]))
    if(any(is.na(isNumeric)) && mutModel=='Stepwise'){
      f_errorWindow("Stepwise mutation model require only numeric allele names. 
                    Change allele names or choose a different mutation model.")
    }
    
    # db2 <- numeric()
    # for(i in 1:ncol(freqs)){
    #   ix <- which(!is.na(freqs[,i]))
    #   a <- alleleNames[ix]
    #   f <- freqs[ix,i]
    #   Marker <- rep(markerNames[i],length(ix))
    #   db2New <- data.frame(Marker,a,f)
    #   colnames(db2New) <- c("Marker","Allele","Frequency")
    #   db2 <- rbind(db2,db2New)
    # }
    #Make lists of alleles and frequencies per marker
    allelesAll <- split(db2$Allele,db2$Marker)
    afreqAll <- split(db2$Frequency,db2$Marker)

    #Pedigree data
    ped1 <- get("ped1",envir=mmTK)
    ped2 <- get("ped2",envir=mmTK)
    persons1 <- get("persons_ped1",envir=mmTK)
    persons2 <- get("persons_ped2",envir=mmTK)
    pedigrees <- list(ped1,ped2)
    idxC <- get("idxC",envir=mmTK) #Individuals who are contributors
    idxK <- unique(G[,1]) #Individuals with known genotypes
    idxU <- idxC[!idxC%in%idxK] #Contributors with uknown genotypes

    #Dropout/drop-in
    #Set default dropout 0 for all contributors if
    if(!exists('dropliste',envir=mmTK)){
      dDef <- get('drop',envir=mmTK)
      D <- rep(list(dDef$d),length(idxC))
      names(D) <- idxC
      di <- dDef$di
    } else{
      drop <- get('dropliste',mmTK)
      D <- drop[idxC]
      di <- drop[["dropin"]]
    }
    

    ############# Computations ###########
    #infoWindowLR <- gwindow("",visible=TRUE)
    #glabel("Computing LR...",container=infoWindowLR)
    
    markers <- names(R)
    LRmarker <- numeric(length(markers))
    for(i in 1:length(R)){
      #Create locus object for each marker
      alleles <- as.character(allelesAll[[which(names(allelesAll)==markers[i])]])
      afreq <- afreqAll[[which(names(afreqAll)==markers[i])]]
      locus <- FamiliasLocus(frequencies=afreq,name=markers[i],
                             allelenames=alleles, MutationModel=mutModel, femaleMutationRate=r1,
                             maleMutationRate=r2, MutationRange=range1)
      
      #If there is a silent allele we need to modify the mutation matrix.
      #Silent allele must be given as 's' (not 'silent' as in Familias)
      #That way Familias will treat it like a regular allele, 
      #while relMix will treat is specially
      if('Silent'%in%alleles){
        newAlleles <- c(alleles[-length(alleles)],'s')
        mm <- locus$femaleMutationMatrix #Assuming same mutation matrix for male and female
        colnames(mm) <- rownames(mm) <- newAlleles
        locus <- FamiliasLocus(frequencies=afreq,name=markers[i],  
                               allelenames= newAlleles, MutationModel='Custom', MutationMatrix=mm)
      }
      
      datamatrix <- createDatamatrix(locus,knownGenos[[which(names(knownGenos)==markers[i])]],idsU=idxU)
      names(pedigrees) <- c("Paternity","Unrelated")
      res <- relMix(pedigrees, locus, R=R[[i]], datamatrix, ids=idxC, D=lapply(D,function(x) c(x,x^2)),di=di, kinship=theta)
      LRmarker[i] <- res$Paternity/res$Unrelated
    }
    
    #dispose(infoWindowLR)
    
    #Data <- data.frame(Marker=c(markers,"Total"),LR=c(LRmarker,prod(LRmarker)))
    Data <- data.frame(Marker=markers,LR=round(LRmarker,4))
    LRwindow <- gwindow("Results",visible=FALSE)
    LRgroup1 <- ggroup(horizontal=FALSE,container=LRwindow)
    glabel(paste("Total LR:",format(prod(LRmarker),,scientific=TRUE)),container=LRgroup1)
    LRgroup2 <- ggroup(container=LRwindow)
    tab <- gtable(Data,container=LRgroup2)
    LRbut <- gbutton("Save LR to file", container=LRwindow,handler=function(h,...){
      f_export(rbind(Data,data.frame(Marker="Total",LR=format(prod(LRmarker),scientific=TRUE))))
      dispose(h$obj)
    })
    visible(LRwindow) <- TRUE
  }


  ############################################################

  mmTK = new.env() #create new environment object

  #Assign some default parameters
  defaultList <- list(Mut='Equal',
                    mutPar=list(mutRange=0.5,
                              fMutRate=0,
                              mMutRate=0),
                    optPar=list(theta=0,
                                silent=0,
                                MAF=0.001),
                    drop=list(d=0,
                              di=0)
                    )
  
  for(i in 1:length(defaultList)){
    assign(names(defaultList)[i],defaultList[[i]],envir=mmTK)
  }
 

  win <- gwindow("RelMix",height=500,width=500,visible=FALSE)
  group1 <- ggroup(horizontal = TRUE, container=win,spacing=10)

  ###### Import data #####
  #Buttons for importing files
  dataGroup <- ggroup(horizontal = FALSE, container=group1,spacing=7)
  impFrame <- gframe("Import data",container=dataGroup,horizontal = FALSE)
  objMix <- gbutton(text="Import mixture profile",container=impFrame,handler=f_importprof,action="mixture")
  objRef <- gbutton(text="Import reference profiles",container=impFrame,handler=f_importprof,action="reference")
  objRef <- gbutton(text="Database",container=impFrame,handler=f_database)

  ###### Mutations #####
  #Mutation frame
  mutFrame <- gframe("Mutations",container=group1,horizontal=FALSE,spacing=7)
  mutGroup <- ggroup(horizontal = FALSE, container=mutFrame,spacing=7)
  mutButton <- gbutton("Mutations", handler=f_mutations, container=mutGroup)

  group2 <- ggroup(horizontal = TRUE, container=win,spacing=10)

  ####### Pedigrees ######
  #Pedigree frame
  pedGroup <- ggroup(horizontal = FALSE, container=group2,spacing=7)
  pedFrame <- gframe("Pedigrees",container=pedGroup)
  #Should we have the option to include different contributors under each hypothesis???
  #Pedigree 1 button
  glabel("Pedigree 1",container=pedFrame)
  objPed <- gcombobox(c("Paternity","Unrelated","Custom"), selected=0,container=pedFrame,handler=f_pedigree,action="1")
  #Pedigree 2 button
  glabel("Pedigree 2",container=pedFrame)
  objPed2 <- gcombobox(c("Paternity","Unrelated","Custom"), selected=0,container=pedFrame,handler=f_pedigree,action="2")

  ######## Contributors #######
  contFrame <- gframe("Contributors",container=group2,horizontal=TRUE)
  objCont <- gbutton(text="Specify contributors",container=contFrame,handler=f_contributors)

  ######## Dropout/drop-in #########
  group3 <- ggroup(horizontal = TRUE, container=win,spacing=10)
  #Dropout button
  dropFrame <- gframe("Dropout and drop-in",container=group3,horizontal=TRUE)
  dropButton <- gbutton("Specify dropout and drop-in", handler=f_dropout, container=dropFrame)

  ####### Compute LR ######
  group4 <- ggroup(horizontal = TRUE, container=win,spacing=10)
  #Compute LR button
  LRbutton <- gbutton("Compute LR", container=group4, handler=function(h,...){
    f_LR()
    })

  visible(win) <- TRUE
}



