
# convert factors.to.strings.in.dataframe 
 convert.factors.to.strings.in.dataframe <- function(dataframe)
    {
        class.data  <- sapply(dataframe, class)
        factor.vars <- class.data[class.data == "factor"]
        for (colname in names(factor.vars))
        {
            dataframe[,colname] <- as.character(dataframe[,colname])
        }
        return (dataframe)
    }
    



