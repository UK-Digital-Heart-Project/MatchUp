# Tim's code for matching three or more patient groups on the basis of covariates
# March 2018, V 1.1, tim.dawes@imperial.ac.uk


# Demo version: HOW TO RUN
#
# 1. Type 'getwd()' [Enter] in the Console to find which folder you are in now
# 2. Put 'Matching.csv' in this folder
# 3. Run the function below by typing: output<- MatchUp("Matching.csv", final.group.size = ...) [Enter]*
# 4. The output is now stored in the variable 'output' and can be saved by typing: write.csv(output, file = "Output.csv") [Enter]
# 5. The output data is the same as the input data with one extra column which denotes whether included in the final matching or not




    


MatchUp<- function (datafile, final.group.size)  {
          options(warn=-1)
          covariates<- read.csv(datafile, header=TRUE)
          ID.column.number<- which(grepl("ID", colnames(covariates))==TRUE)
          
        # Error messages
          if (length(which(colnames(covariates)==toupper(colnames(covariates))))==0) {return("\n Covariates in uploaded datafile must be indicated by putting the column in upper case please.")}
          if (sum(ID.column.number)==0) {return("Please include the letters 'ID' in exactly one column to mark the subject identifiers")}
            
        # Scale the data so each covariate is equally weighted in OLS
          covariates.ORIGINAL<- covariates
          match.on.these.cols<- which(colnames(covariates)==toupper(colnames(covariates)))
          sc<- ce<- rep(0, length(match.on.these.cols))
          for (i in match.on.these.cols) {ce[i]<- mean(covariates[,i]); sc[i]<- sd(covariates[,i]); covariates[,i]<- (covariates[,i]-ce[i])/sc[i]}
        
        # Add a group number to each group and define colours for the plotting
          numbergroups<- length(unique(covariates$Group))
    
          
        # Define variables
          diff<- list()
          means<- matrix(0, nrow=numbergroups, 4)
          hit.minimum.group.size<- rep(FALSE,numbergroups)
          iter<- 1
          covariates$NEWGROUP<- rep(0,nrow(covariates))
          
          
          repeat
            {
              for (target.group in which(hit.minimum.group.size==FALSE))
                {
                  pooled.groups<- which(hit.minimum.group.size==FALSE)[-target.group]
                  
                  covariates$NEWGROUP<- rep(0, nrow(covariates))
                  covariates$NEWGROUP[which(covariates$Group %in% target.group == TRUE)]<- 1
                  
                  if (sum(covariates$NEWGROUP)<=final.group.size) {hit.minimum.group.size[target.group]<- TRUE; break}
                  if (sum(hit.minimum.group.size[-target.group])==length(hit.minimum.group.size[-target.group])) {no2takeaway<- length(which(covariates$Group==target.group))-final.group.size} else {no2takeaway=1}
                  
                  f<- paste("NEWGROUP", paste(colnames(covariates)[match.on.these.cols], collapse=" + "), sep=" ~ ")
                  m <- glm(as.formula(f), data = covariates, family=binomial())
                  prs<- data.frame(ps = predict(m, type = "response"), NEWGROUP = m$model$NEWGROUP)
                  
                  o<- order(prs$ps[which(prs$NEWGROUP==1)],decreasing=TRUE)
                  t<- covariates[which(prs$NEWGROUP==0),]
                  s<- covariates[which(prs$NEWGROUP==1)[o[-(1:no2takeaway)]],]
                  covariates<- rbind(t,s)
                
                }
            if (sum(hit.minimum.group.size)==length(hit.minimum.group.size)) {break}
            iter<- iter + 1
            }
          
          covariates.ORIGINAL$MATCHEDGROUP<- !is.na(match(covariates.ORIGINAL[,ID.column.number], covariates[,ID.column.number]))
          return (covariates.ORIGINAL)
      }
      
      
      

      

# *Notes:
# Matching.csv must be (a) .csv format, (b) column names for each column, (c) all covariate names in capitals / no non-covariates in capitals, (d) ID columnname must contain the letters, in capitals, "ID", for example:
# "Unique ID", "Unique.ID", "IDs" would all be fine; "uniqueid","uniqueId","ids" would not be fine.
# Final.group.size decides the number in each final group. Small samples may lead to non-convergence in the glm, although the warnings for this have been switched off


# Assumptions:
# 1. There are three or more non-overlapping clinical groups which need to be matched on the basis of covariates
# 2. The initial groups may be of unequal sizes
# 3. You want to match the groups 1:1:1
# 4. There are no missing data


# Design:
# i. Propensity matching approach
# ii. Sequential comparison of each group ('target' group) with the remaining groups ('pooled' groups)
# iii. Subject in the target group with lowest propensity score is discarded, provided the number of subjects in the target group exceeds a minimum threshold


