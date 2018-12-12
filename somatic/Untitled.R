
sample <- read.table("/Users/gdemidov/Conferences_and_workshops/2018/MSc_course/Exercises/Assignment_04_Rosanna_Krebs/sample2.2.txt")
population_params <- read.table("/Users/gdemidov/Conferences_and_workshops/2018/MSc_course/Exercises/Assignment_04_Rosanna_Krebs/population_parameter2.1.txt")

sampleZ <- (sqrt(sample) - population_params[,1]) / population_params[,2]
plot(sampleZ[,1])


myIris <- iris
myIris$smallSepalLength <- round(myIris$Sepal.Length/5,2)
#myIris$smallSepalLength <- myIris$Sepal.Length/5 # wrong answer
mean(myIris[!myIris$smallSepalLength %in% c(1.08,0.98,0.94,0.92,1.00),"Sepal.Width"])
tol = 10**-10
mean(myIris[!all( abs(myIris$smallSepalLength - c(1.08,0.98,0.94,0.92,1.00)) < tol ) ,"Sepal.Width"])
