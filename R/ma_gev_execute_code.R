# Download and install the magev package from the attached magev_0.0.0.9000.tar file
# Then, load the package using library(magev)
library(magev)  

# Load the Haenam dataset  
data(haenam)  

# Apply the ma.gev.open function
ma.gev.open(data=haenam[,-1], quant=c(.98, .99, .995), weight='like0')
