if(Ds|Dr){
  if(Ds&Dr){
    coln <- colnames(DrInput_Target)[-9]
    Target1 <- unique(merge(DisInput_Target,DrInput_Target, by=coln,all=T))
  }else{
    if(Ds)Target1 <- DisInput_Target
    if(Dr)Target1 <- DrInput_Target
  }
  if(length(Target1$ORGANISM.x)>0)Target1$ORGANISM.x <- factor(Target1$ORGANISM.x)
  DesiredAct <- unique(Target1[!is.na(Target1$Action.x),c('GENE','Action.x')]) #16
  colnames(DesiredAct) <- c('GENE','DesiredAct')
  if(length(DesiredAct$GENE)>0)Target1 <- merge(Target1,DesiredAct,by='GENE',all.x=T)
  
} #276

if(Gn|Cs){
  if(Gn&Cs){
    Target2 <- unique(merge(Custom_Target,GeneInput_Target,by='GENE',all=T))
  }else{
    if(Gn)Target2 <- GeneInput_Target
    if(Cs)Target2 <- Custom_Target
  }
} #25

T1 <- length(Target1$GENE)
T2 <- length(Target2$GENE)

if(T1|T2){
  if(T1&T2){
    Result1_Target <- unique(merge(Target1,Target2[!is.na(Target2$GENE),],by='GENE',all=T))
  }else{
    if(T1)Result1_Target <- Target1
    if(T2)Result1_Target <- Target2
  }
} #276


#### Combined Dataset Procesing

# Find Drug Candidate
Result2_DrugCand <- unique(merge(Result1_Target,DBDC_Target_Full,by='GENE',all.x=T)) #8737
Result2_DrugCand$Disease <- input$checkbox_DiseaseName
Result2_DrugCand$Drug_Found <- FALSE
Result2_DrugCand[!is.na(Result2_DrugCand$Drug),'Drug_Found'] <- TRUE
Result2_DrugCand$ORGANISM <- factor(Result2_DrugCand$ORGANISM)

# Find Known Indication of Drug Candidates
Result3_Ix <- merge(Result2_DrugCand,DrIx_byDrug,by.x='Drug',by.y='DRUG_NAME',all.x=T) #8737

# Find Target-associated Disease(Condition)
#Result3_Ix <- merge(Result3_Ix,Dis_byTarget,by='GENE',all.x=T) #8737

# Remove Duplicated (Indicated Drug & Candidate Drug) Rows
# Assess Same Indication, Same Action of Drug on Target
Result4 <- Result3_Ix
if(nrow(Result4[!is.na(Result4$Drug.x),])>0){
  Result4 <- unique(Result4[!(!is.na(Result4$Drug.x)&Result4$Drug.x==Result4$Drug),]) #8439
  if(length(Result4$Disease)>0) Result4 <- Result4[!is.na(Result4$Disease),] #18866
  Result4$Same_Act <- apply(Result4[,c('Action.x','Action')],1,function(x){x[1]==x[2]})
  if(length(Result4$DesiredAct)>0){
    Result4$Same_Act_expected <- apply(Result4[,c('DesiredAct','Action')],1,function(x){x[1]==x[2]})
  }
  Result4$Same_Organism <- apply(Result4[,c('ORGANISM.x','ORGANISM')],1,function(x){x[1]==x[2]})
  Result4$Same_Indication <- FALSE
  Result4[Result4$Drug %in% Result4$Drug.x,'Same_Indication'] <- TRUE
}
#Result4$Same_Act <- apply(Result4[,c('Action.x','Action.y')],1,function(x){x[1]==x[2]})

# Match Drug Group (Approval status)
Result4 <- merge(Result4,Drug_Group,by='DBID',all.x=T)

ColnameSet <- data.frame(SelectCol=c('Disease','Drug.x','GENE','Drug','Known_Indications','Targ_Frq_Dr',
                                     'Same_Act','Same_Act_expected','Same_Organism','Same_Indication','Action.x','Action','ORGANISM.x','ORGANISM',
                                     'Drug_Found','DBID','Source_DC.x','Source_DB.x','Source_DGN','Source_Custom',
                                     'Source_DC','Source_DB','Input_Disease','Input_Drug','Input_Gene',
                                     'Pharm_active','Approved','Investigational','Experimental','Withdrawn','Nutraceutical','Small','Biotech','Illicit'),
                         RenameCol=c('Disease','Drug_Indicated','Target_Gene','Drug_Candidate','Known_Indications','Targ_Frq_Dr',
                                     'Same_Act','Same_Act_expected','Same_Organism','Same_Indication','Action_x','Action_y','Species_x','Species_y',
                                     'Drug_Found','DBID','Source_DC','Source_DB','Source_DGN','Source_Custom',
                                     'Source_Drug_DC','Source_Drug_DB','Input_Diseaes','Input_Drug','Input_Gene',
                                     'Pharm_active','Approved','Investigational','Experimental','Withdrawn','Nutraceutical','Small','Biotech','Illicit'))

Coldiff <- setdiff(ColnameSet$SelectCol,colnames(Result4))
Col_select <- setdiff(ColnameSet$SelectCol,Coldiff)

Result5 <- unique(Result4[,Col_select]) #7178

for(cl in Col_select[grepl('Source_|Input_',Col_select,ignore.case = T)]){
  Result5[is.na(Result5[,cl]),cl] <- FALSE
}

colnames(Result5) <-  ColnameSet[ColnameSet$SelectCol %in% Col_select,'RenameCol']

return(Result5)