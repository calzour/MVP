"""
error calculation
takes in a list and finds the uncertainty from variances
requires VarianceOfList() function
"""

def bootstrapErrors(List, N, b):
    #calculate parameters
    noSamples = 100
    listSize = len(List)
    sampleList = [] #empty list for sampling
    obsList = [] #list for variances of sample lists
    variance = 0.0 #the error that the function returns
    #dummy variable for variance function
    dummy = 0.0

    for p in range(noSamples): #create 100 random samples
        sampleList = np.random.choice( List, listSize )
        tempVariance, dummy = varianceOfList(sampleList)

        obsList.append(tempVariance)
    #after 100 times, find the variance of this new list
    varianceSq, dummy = varianceOfList(obsList)

    variance = np.sqrt(varianceSq)

    return(variance)
